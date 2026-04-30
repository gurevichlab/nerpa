use std::fs;
use std::path::{Path, PathBuf};
use anyhow::{anyhow, Result};

#[derive(Debug)]
struct SvgInfo {
    // Logical coordinate size we’ll use for scaling/stacking
    width: f32,
    height: f32,

    // viewBox offset (often 0,0). If present, we’ll translate to normalize.
    min_x: f32,
    min_y: f32,

    // Everything inside the outer <svg>...</svg>
    inner: String,
}

fn parse_leading_f32(s: &str) -> Option<f32> {
    // Accept things like "100", "100px", "100.5pt"
    let mut end = 0usize;
    for (i, ch) in s.char_indices() {
        if ch.is_ascii_digit() || ch == '.' || ch == '-' || ch == '+' {
            end = i + ch.len_utf8();
        } else {
            break;
        }
    }
    if end == 0 {
        return None;
    }
    s[..end].parse::<f32>().ok()
}

fn parse_attrs(opening_tag: &str) -> std::collections::HashMap<String, String> {
    // Tiny attribute parser for key="value" pairs.
    // Assumes attributes are quoted with "..." (common for SVG).
    let mut attrs = std::collections::HashMap::new();
    let bytes = opening_tag.as_bytes();
    let mut i = 0usize;

    // Skip until after "<svg"
    while i < bytes.len() && bytes[i] != b' ' && bytes[i] != b'\t' && bytes[i] != b'\n' {
        i += 1;
    }

    while i < bytes.len() {
        // Skip whitespace and tag end chars
        while i < bytes.len()
            && (bytes[i] == b' ' || bytes[i] == b'\t' || bytes[i] == b'\n' || bytes[i] == b'>' || bytes[i] == b'/')
        {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }

        // Read key
        let key_start = i;
        while i < bytes.len() && bytes[i] != b'=' && bytes[i] != b' ' && bytes[i] != b'\t' && bytes[i] != b'\n' && bytes[i] != b'>' {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'=' {
            // Not a key=value; bail out
            break;
        }
        let key = opening_tag[key_start..i].trim().to_string();
        i += 1; // '='

        // Expect quote
        if i >= bytes.len() || bytes[i] != b'"' {
            // Only handling double-quoted attrs
            break;
        }
        i += 1; // opening quote

        let val_start = i;
        while i < bytes.len() && bytes[i] != b'"' {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        let val = opening_tag[val_start..i].to_string();
        i += 1; // closing quote

        if !key.is_empty() {
            attrs.insert(key, val);
        }
    }

    attrs
}

fn parse_svg_info(svg_text: &str) -> Result<SvgInfo> {
    let svg_start = svg_text
        .find("<svg")
        .ok_or_else(|| anyhow!("No <svg ...> start tag found"))?;

    let tag_end = svg_text[svg_start..]
        .find('>')
        .map(|idx| svg_start + idx)
        .ok_or_else(|| anyhow!("No closing '>' for <svg ...> tag"))?;

    let opening_tag = &svg_text[svg_start..=tag_end];

    let close_tag = "</svg>";
    let svg_end = svg_text
        .rfind(close_tag)
        .ok_or_else(|| anyhow!("No </svg> end tag found"))?;

    let inner = svg_text[tag_end + 1..svg_end].to_string();

    let attrs = parse_attrs(opening_tag);

    // Prefer viewBox for logical sizing when available (more reliable than width="100%").
    let mut min_x = 0.0f32;
    let mut min_y = 0.0f32;
    let mut vb_w: Option<f32> = None;
    let mut vb_h: Option<f32> = None;

    if let Some(view_box) = attrs.get("viewBox") {
        let parts: Vec<&str> = view_box.split_whitespace().collect();
        if parts.len() == 4 {
            min_x = parts[0].parse::<f32>().unwrap_or(0.0);
            min_y = parts[1].parse::<f32>().unwrap_or(0.0);
            vb_w = parts[2].parse::<f32>().ok();
            vb_h = parts[3].parse::<f32>().ok();
        }
    }

    let width = if let (Some(w), Some(_h)) = (vb_w, vb_h) {
        w
    } else {
        attrs
            .get("width")
            .and_then(|s| parse_leading_f32(s))
            .ok_or_else(|| anyhow!("Missing/invalid width (and no usable viewBox)"))?
    };

    let height = if let (Some(_w), Some(h)) = (vb_w, vb_h) {
        h
    } else {
        attrs
            .get("height")
            .and_then(|s| parse_leading_f32(s))
            .ok_or_else(|| anyhow!("Missing/invalid height (and no usable viewBox)"))?
    };

    Ok(SvgInfo {
        width,
        height,
        min_x,
        min_y,
        inner,
    })
}


/// Stacks multiple SVG files vertically (one on top of another), scaling each
/// so they all match the widest SVG's logical width, and writes the result.
///
/// - `input_paths`: SVGs to stack, in order (top to bottom)
/// - `output_path`: where to write the combined SVG
pub fn join_svgs_vertical(input_paths: &[&Path], output_path: &Path) -> Result<()> {
    if input_paths.is_empty() {
        return Err(anyhow!("input_paths is empty"));
    }

    let mut svgs = Vec::with_capacity(input_paths.len());
    for p in input_paths {
        let txt = fs::read_to_string(p)?;
        let info = parse_svg_info(&txt)?;
        svgs.push(info);
    }

    let max_width = svgs
        .iter()
        .map(|s| s.width)
        .fold(0.0f32, |a, b| a.max(b));

    if max_width <= 0.0 {
        return Err(anyhow!("Computed max width is <= 0"));
    }

    // Compute total height after scaling each to max_width
    let mut total_height = 0.0f32;
    let mut placements: Vec<(f32, f32)> = Vec::with_capacity(svgs.len()); // (y_offset, scale)

    for s in &svgs {
        let scale = max_width / s.width;
        let scaled_h = s.height * scale;
        placements.push((total_height, scale));
        total_height += scaled_h;
    }

    let mut out = String::new();
    out.push_str(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    out.push('\n');
    out.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"#,
        w = max_width,
        h = total_height
    ));
    out.push('\n');

    for (s, (y, scale)) in svgs.iter().zip(placements.iter()) {
        // Transform order matters (applied right-to-left):
        // 1) translate(-min_x, -min_y) to normalize viewBox origin (if any)
        // 2) scale to match max width
        // 3) translate down by y to stack
        out.push_str(&format!(
            r#"<g transform="translate(0 {y}) scale({scale}) translate({tx} {ty})">"#,
            y = y,
            scale = scale,
            tx = -s.min_x,
            ty = -s.min_y,
        ));
        out.push('\n');
        out.push_str(&s.inner);
        out.push('\n');
        out.push_str("</g>\n");
    }

    out.push_str("</svg>\n");

    fs::write(output_path, out)?;
    Ok(())
}

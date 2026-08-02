pub struct SplitWhen<I, F>
where
    I: Iterator,
{
    iter: I,
    pred: F,
    pending: Option<I::Item>,
}

pub fn split_when<I, F>(iterable: I, pred: F) -> SplitWhen<I::IntoIter, F>
where
    I: IntoIterator,
    F: FnMut(&I::Item, &I::Item) -> bool,
{
    SplitWhen {
        iter: iterable.into_iter(),
        pred,
        pending: None,
    }
}

impl<I, F> Iterator for SplitWhen<I, F>
where
    I: Iterator,
    F: FnMut(&I::Item, &I::Item) -> bool,
{
    type Item = Vec<I::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        let first = match self.pending.take().or_else(|| self.iter.next()) {
            Some(x) => x,
            None => return None,
        };

        let mut chunk = vec![first];

        while let Some(item) = self.iter.next() {
            let prev = chunk.last().unwrap();
            if (self.pred)(prev, &item) {
                self.pending = Some(item); // start next chunk with this item
                break;
            }
            chunk.push(item);
        }

        Some(chunk)
    }
}

pub fn rounded(f: f64, digits: usize) -> f64 {
    if !f.is_finite() {
	return f;
    }
    if digits == 0 {
	return f.round();
    }
    let factor = 10f64.powf(digits as f64);
    if !factor.is_finite() || factor == 0.0 {
	return f;
    }
    (f * factor).round() / factor
}

use std::{borrow::Borrow, cmp::Ordering};

pub fn f64_cmp<T: Borrow<f64>>(a: &T, b: &T) -> Ordering {
    let f1 = a.borrow();
    let f2 = b.borrow();

    if f1.is_nan() && f2.is_nan() {
	std::cmp::Ordering::Equal
    } else if f1.is_nan() {
	std::cmp::Ordering::Greater
    } else if f2.is_nan() {
	std::cmp::Ordering::Less
    } else {
	f1.partial_cmp(f2).unwrap()
    }
}

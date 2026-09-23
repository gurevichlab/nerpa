// -- spectrum --
let spectrum;
let spectrumCanvas;

function drawSpectrum(spectrumID, variantID){
    document.getElementById('spectrumCanvas').innerHTML = "";

    const rect = document.getElementById('spectrum').getBoundingClientRect()
    const canvasWidth = Math.round(rect.width * 0.9);
    const canvasHeight = Math.round(rect.height * 0.5);
    spectrumCanvas = new ChemDoodle.PerspectiveCanvas('spectrumCanvas', canvasWidth, canvasHeight);
    spectrumCanvas.styles.plots_color="grey";
    spectrumCanvas.styles.plots_width= 1;
    spectrumCanvas.styles.text_font_size = 14;
    spectrumCanvas.styles.text_font_families = ["Arial", "Charcoal", "sans-serif"];

    const peakArray = spectra[spectrumID].peaks;
    const spectrumJcampFile = peakArrayToJcamp(peakArray);
    spectrum = ChemDoodle.readJCAMP(spectrumJcampFile); 
    spectrumCanvas.loadSpectrum(spectrum); 
    colorMatchedPeaks(spectrumID, variantID);

    const oldRepaint = spectrumCanvas.repaint;
    spectrumCanvas.repaint = function(e) {
        oldRepaint.call(this,e);
        updatePeaks();
    };  

    const svgSpectrum = document.getElementById('spectrumSvg');
    svgSpectrum.childNodes.forEach(peakLine => {
        peakLine.addEventListener('click', e => {
            const nid = peakLine.getAttribute('nid').split(',');
            // Important, because peaks are highlighted differently depending on whether they were selected directly from a peak/peak row or from a monomer.
            peakLine.setAttribute("clicked", true);
            selectMS(nid);
        });
    });
    
    const infoDiv = document.getElementById("spectrumInfo");
    const match = getSpectrumVariantMatch(spectrumID, variantID);

    let infoText = `<span style="color: #000000; line-height: 2;"><strong>Spectrum ID: ${spectrumID}</strong></span>
                    <br> <span><strong>Precursor:</strong> ${spectra[spectrumID].precursor_mass.toFixed(3)} Da</span>
                    <span style="margin: 0 10px; color: #bbb;">|</span> 
                    <span><strong>MS Score:</strong> ${match.score}</span>`;
            
    if(spectra[spectrumID].charge != null){
        infoText += `<span style="margin: 0 10px; color: #bbb;">|</span> 
                     <span><strong>Charge:</strong> +${spectra[spectrumID].charge}</span>`;
    } 
    infoDiv.innerHTML = infoText;

    initPeakTable(svgSpectrum);
    window.addEventListener("resize", (event) => {
        resizeSpectrum();
    }); 
}
function resizeSpectrum(){
    if(spectrumCanvas === undefined) return;
    const rect = document.getElementById('spectrum').getBoundingClientRect()
    const canvasWidth = Math.round(rect.width * 0.9);
    const canvasHeight = Math.round(rect.height * 0.5);
    spectrumCanvas.resize(canvasWidth,canvasHeight);
    spectrumCanvas.loadSpectrum(spectrum); 
}
function initPeakTable(svgSpectrum){
    const tbPeak = document.querySelector("#peakTable tbody");
    
    tbPeak.innerHTML = '';
    for(const peak of svgSpectrum.childNodes){
        if(!peak.id) continue;

        const line = tbPeak.insertRow();
        line.setAttribute("nid", peak.getAttribute("nid"));
        line.addEventListener('click', e => {
            const nid = peak.getAttribute("nid").split(',');
            // Important, because peaks are highlighted differently depending on whether they were selected directly from a peak/peak row or from a monomer.
            peak.setAttribute("clicked", true);
            selectMS(nid);
        });
        const peakMzCell = line.insertCell();
        peakMzCell.innerHTML = peak.getAttribute('mz');  

        const peakIntensityCell = line.insertCell();
        peakIntensityCell.innerHTML = peak.getAttribute('intensity');
        
        const chargeCell = line.insertCell();
        chargeCell.innerHTML = peak.getAttribute('charge');  
  
        const peakMassErrorAbsCell = line.insertCell();
        peakMassErrorAbsCell.innerHTML = peak.getAttribute('massErrorAbs');  

        const peakMassErrorRelCell = line.insertCell();
        peakMassErrorRelCell.innerHTML = peak.getAttribute('massErrorRel');
    }
}
function peakArrayToJcamp(peaks){
    let jcamp = `##TITLE=${''}\n##DATA TYPE= MASS SPECTRUM\n##XUNITS= m/z\n##YUNITS= relative abundance\n##PEAK TABLE= (XY..XY)\n`;
    for (const p of peaks) {
        const x = p.mass;
        const y = p.intensity;
        jcamp = jcamp + `${x}, ${y}\n`
    }
    return jcamp + '##END=\n';
}
function colorMatchedPeaks(spectrumID, variantID){
    const svgSpectrum = document.getElementById('spectrumSvg');
    svgSpectrum.innerHTML = "";

    const variantObject = getVariantObject(variantID);
    const mod_map = variantObject.old_to_new_mon_map;
    const match = getSpectrumVariantMatch(spectrumID, variantID);

    const specWidth = spectrum.memory.width;
    const specHeight = spectrum.memory.height;
    const specOffsetLeft = spectrum.memory.offsetLeft;
    const specOffsetBottom = spectrum.memory.offsetBottom;
    const specOffsetTop = spectrum.memory.offsetTop;

    const origY = spectrum.getTransformedY(0,  spectrumCanvas.styles, specHeight, specOffsetBottom, specOffsetTop);
  
    for(const peak of spectrum.data){
        const matchedPeak = getMatchObject(match.matched_peaks, peak.x);
        if(matchedPeak === undefined) continue; // i.e. peak is not matched 
      
        const coordx = spectrum.getTransformedX(peak.x, spectrumCanvas.styles, specWidth,specOffsetLeft);
        const coordy = spectrum.getTransformedY(peak.y,  spectrumCanvas.styles, specHeight, specOffsetBottom, specOffsetTop);
      
        const nid = translateMask(matchedPeak.theoretical_fragment_mask, mod_map);

        const mz = peak.x;
        const intensity = spectra[spectrumID].peaks[matchedPeak.experimental_peak_idx].intensity; //peak.y is abs abundance not intensity
        const charge = matchedPeak.charge;
        const massErrorAbs = parseFloat(matchedPeak.theoretical_mz - matchedPeak.experimental_mz);
        const massErrorRel = parseFloat((massErrorAbs / matchedPeak.experimental_mz ) * 1000000);
 
        const peakLine = document.createElementNS(svgns, 'line');
        peakLine.setAttribute("id", `${peak.x}`);
        peakLine.setAttribute("nid", `${nid}`);
        peakLine.setAttribute("clicked", false);
        // set peak info to prepare for initialization of peak table.
        peakLine.setAttribute("mz", `${mz.toFixed(3)}`);
        peakLine.setAttribute("intensity", `${intensity.toFixed(1)}`);
        peakLine.setAttribute("charge", `${charge}`);
        peakLine.setAttribute("massErrorAbs", `${massErrorAbs.toFixed(3)}`);
        peakLine.setAttribute("massErrorRel", `${massErrorRel.toFixed(1)}`);

        peakLine.setAttribute('x1', coordx);
        peakLine.setAttribute('y1', origY);
        peakLine.setAttribute('x2', coordx);
        peakLine.setAttribute('y2', coordy);
        peakLine.setAttribute('stroke', '#2a6881');
        peakLine.setAttribute('stroke-width', '4');
        svgSpectrum.appendChild(peakLine);
    }
}
function updatePeaks(){
    
    const specWidth = spectrum.memory.width;
    const specHeight = spectrum.memory.height;
    const specOffsetLeft = spectrum.memory.offsetLeft;
    const specOffsetBottom = spectrum.memory.offsetBottom;
    const specOffsetTop = spectrum.memory.offsetTop;

    const origY = spectrum.getTransformedY(0,  spectrumCanvas.styles, specHeight,specOffsetBottom, specOffsetTop);
    const origX = spectrum.getTransformedX(0,  spectrumCanvas.styles, specWidth, specOffsetLeft);

    for(const peak of spectrum.data){

        const coordx = spectrum.getTransformedX(peak.x, spectrumCanvas.styles, specWidth, specOffsetLeft);
        const coordy = spectrum.getTransformedY(peak.y,  spectrumCanvas.styles, specHeight, specOffsetBottom, specOffsetTop);

        const peakLine = document.getElementById(`${peak.x}`);
        if(peakLine === null || coordx < origX) continue;       // peak is not matched or outside of visible area 

        peakLine.setAttribute('x1', coordx);
        peakLine.setAttribute('y1', origY);
        peakLine.setAttribute('x2', coordx);
        peakLine.setAttribute('y2', coordy);
    }
}
function translateMask(mask, mod_map){
    const deletedIDs = mod_map.filter(t => t[0] != null && t[1] === null).map(t => t[0]);
    const results = [];
    let id = 1;
    for(const char of mask.split('')){
        if(deletedIDs.includes(id)){
            id++;
        }
        const digit = parseInt(char);
        if(digit === 1){
            results.push(id);
        }
        id++;
    }

    return results;
}
function getSpectrumVariantMatch(spectrumID, variantID){
    return spectra_matching_results.find(
            (entry) => 
                entry.spectrum_id === spectrumID &&
                entry.structure_id === variantID
            );
}
function getMatchObject(matchedPeaks, mz){
    return matchedPeaks.find(p => p.experimental_mz === mz);   
}
// -- Nerpa MS modification graph --
function drawOrigGraph(nrpID, variant){
   
    drawGraph(nrpID);

    const variantObj = getVariantObject(variant);
    const mod_map = variantObj.old_to_new_mon_map;
    network.nodes.forEach(n => network.nodes.update(
        { id: n.id, color: getNodeColorOrig(n.id, mod_map) }
    ));
}

function drawModGraph(nrpID, variantID) {

    const graphData = monomer_graph_variants[variantID];
  
    const vis_network = buildVisGraph(graphData.nodes, graphData.edges, "graphModNew", nrpID);

    const variantObj = getVariantObject(variantID);
    const mod_map = variantObj.old_to_new_mon_map;
    vis_network.nodes.forEach(n => vis_network.nodes.update(
        { id: n.id, color: getNodeColorNew(n.id, mod_map), fixed: true, chosen: false}
    ));

    vis_network.addEventListener('click',  e => {
        deselect();
        vis_network.fit();
        vis_network.selectNodes([]);  
    });
    
    const p = document.getElementById("graphModP");
    const variantIDs = extractIDs(variantID);
    p.innerHTML =`Modified Monomer Graph <span style="margin: 0 10px; color: #bbb;">|</span> 
                    <span>rank: ${variantIDs.rank}</span>
                    <span style="margin: 0 10px; color: #bbb;">|</span> 
                    <span>number of modifications: ${variantIDs.numMods}</span>
                `; 
}


function getNodeColorOrig(id, mod_map){

    const tup = mod_map.find(tup => tup[0] === parseInt(id));
    if(!tup){
        return "#cbcbcb"
    } else if(tup[1] === null){
        //(old_idx, null) -- the monomer old_idx was removed
        return "#e79898"    
    } else if(tup[0] != null){
        //(old_idx, new_idx) -- the monomer old_idx was substituted (in this case old_idx=new_idx)
        return "#ffff9e"    
    } else {
        return "#cbcbcb"
    }
}

function getNodeColorNew(id, mod_map){
    
    const tup = mod_map.find(tup => tup[1] === parseInt(id));
    if(!tup){
        return "#cbcbcb"
    } else if(tup[0] === null){
        //(null, new_idx) -- the monomer new_idx was inserted
        return "#98e799" 
    } else  if(tup[1] != null){
        //(old_idx, new_idx) -- the monomer old_idx was substituted (in this case old_idx=new_idx)
        return "#ffff9e"    
    } else {
        return "#cbcbcb"
    }
}

function getVariantObject(variantID){
    const variantIDs = extractIDs(variantID);
    const [variantkey, variantObject] = Object.entries(candidate_NRPs).find(([key]) => {
        const entryIDs = extractIDs(key);
        return entryIDs.nrpID === variantIDs.nrpID &&
            entryIDs.bgcID === variantIDs.bgcID;
    });
    return variantObject.new_variants[variantID]
}

// -- variant Graph ---
let variantNetwork;
function drawVariantGraph(nrpID, variantID, spectrumID){
    const graphData = monomer_graph_variants[variantID];
  
    const vis_network = buildVisGraph(graphData.nodes, graphData.edges, "variantGraphImage", nrpID);
    variantNetwork = vis_network;
    
    variantNetwork.addEventListener('click',  e => {
        if(e.nodes.length === 0){ 
            deselect();
            variantNetwork.fit();
        } else {
            selectMS(e.nodes);
        }
    });

    nodeIdLabel.clear();
    for(const entry of graphData.nodes){
        nodeIdLabel.set(entry.id, entry.label);
    }

    displayVarquestMod(spectrumID, variantID);
}

function displayVarquestMod(spectrumID, variantID){
    const match = getSpectrumVariantMatch(spectrumID, variantID);
    const variantObj = getVariantObject(variantID);
    const mod_map = variantObj.old_to_new_mon_map;
    const vidToNid = varquestIdToNid(variantObj.linearization.length, mod_map); //varquest (now kakapo) reindexes indeces of monomer graph 

    let modSum = 0;
    let subText = '';
    for(const mod of match.modifications){
        if(subText != '') subText += ' + ';

        const initMass = mod.initial_monomer_mass;
        const mass_diff = mod.mass_difference;
        modSum += (initMass + mass_diff);

        const modified_node = variantNetwork.nodes.get(vidToNid[mod.monomer_idx].toString());

        let mass_diff_string;
        if (mass_diff < 0){
            mass_diff_string = `- ${mass_diff.toFixed(3) * -1} Da`;
        } else{
            mass_diff_string = `+ ${mass_diff.toFixed(3)} Da`;
        }
        subText += `<span style="
                            background-color: ${modified_node.color};
                            background-clip: padding-box;
                            border: 4px dashed ${darkenColor(modified_node.color)};
                            border-radius: 999px;
                            padding: 2px 6px;
                            ">${modified_node.label}</span>
                            ${mass_diff_string}`;

        variantNetwork.nodes.update({
            id: modified_node.id,
            label: `${modified_node.label}<i>${mass_diff_string}</i>`,
            shapeProperties: {
                borderDashes: [5, 5] 
            },
            borderWidth: 8,  
            borderWidthSelected: 8
        });
    }

    const modLine = document.getElementById('modLine');
    modLine.innerHTML = `${subText} = ${modSum.toFixed(3)} Da`;
}

function darkenColor(hex, percent = 10) {
  hex = hex.replace(/^\s*#|\s*$/g, '');
  
  const num = parseInt(hex, 16);
  let r = (num >> 16) - Math.round(255 * (percent / 100));
  let g = ((num >> 8) & 0x00FF) - Math.round(255 * (percent / 100));
  let b = (num & 0x0000FF) - Math.round(255 * (percent / 100));

  r = Math.max(0, r);
  g = Math.max(0, g);
  b = Math.max(0, b);

  return `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1).toUpperCase()}`;
}

function varquestIdToNid(length, mod_map){
    // varquestID -> nid
    const result = {};
    let nidCounter = 1;
    const deletedIDs = mod_map.filter(t => t[0] != null && t[1] === null).map(t => t[0]);
   
    for (let id = 0; id < length; id++) {
        if(deletedIDs.includes(nidCounter)){
            nidCounter ++;
        }
        result[id] = nidCounter;
        nidCounter ++;
    }
    return result;
}

function selectMS(nid){
    let nidString = nid.map(String);
    let nidInt = nid.map(Number);
    const hasInvalidIdx = nidInt.some(Number.isNaN);


    // select peaks
    const svgSpectrum = document.getElementById('spectrumSvg');
    const spectrumClick =  Array.from(svgSpectrum.childNodes).some(peakLine => peakLine.getAttribute('clicked') === 'true');
    for(const peak of svgSpectrum.childNodes){
        if(!peak.id) continue;
        const nidPeak = peak.getAttribute('nid').split(',');
        const sameNids = nidPeak.length === nidString.length && nidPeak.every(val => nidString.includes(val));
        const containNids = nidPeak.some(val => nidString.includes(val));
        if(sameNids){
            // all nids are explained by this peak 
            peak.setAttribute("opacity", "1");
            peak.removeAttribute('stroke-dasharray');
        } else if (containNids && !spectrumClick){
            // not all nids are explained by this peak and peak was selected via variant graph
            peak.setAttribute("opacity", "1");
            peak.setAttribute('stroke-dasharray', '2 2');
        } else {
            peak.setAttribute("opacity", "0.5");
            peak.removeAttribute('stroke-dasharray');
        }
        // ensures that peak is reset to default when selected again 
        peak.setAttribute('clicked', false);
    }
     
    // select peak Table
    const peakTableBody = document.querySelector('#peakTable tbody');
    [...peakTableBody['rows']].forEach(r => {
        const nidRow = r.getAttribute('nid').split(',');
        const sameNids = nidRow.length === nidString.length && nidRow.every(val => nidString.includes(val));
        if (sameNids){
            r.style.border = "2px solid #9EC37B";
            r.style.background = '#ddfcdb';
            r.scrollIntoView({
                behavior: 'smooth',
                block: 'center'
            });
        } else {
            r.style.cssText = "";
        }
    });
    
    // select variant graph 
    if(hasInvalidIdx) {
        nidString = nidInt.filter(nid => !Number.isNaN(nid)).map(String);
    } 
    if(nidString.length === 0){
        // deselect variant graph 
        variantNetwork.selectNodes([]);  
        increaseTranparency(nidInt, variantNetwork.nodes, variantNetwork.edges);
    } else {
        // select variant graph 
        variantNetwork.selectNodes(nidString);  
        increaseTranparency(nidString, variantNetwork.nodes, variantNetwork.edges);
    } 

    // select variant molecule
    selectMol(nidString);

}
function deselectMS(){

    // deselect variant graph
    variantNetwork.selectNodes([]); 
    increaseTranparency([], variantNetwork.nodes, variantNetwork.edges);

    // deselect peaks
    const svgSpectrum = document.getElementById('spectrumSvg');
    for(const peak of svgSpectrum.childNodes){
        if(!peak.id) continue;
        peak.setAttribute("opacity", "1");
        peak.removeAttribute('stroke-dasharray');
        peak.setAttribute('clicked', false);
    }

    // deselect peak Table
    const peakTable = document.querySelector('#peakTable tbody');
    [...peakTable['rows']].forEach(r => r.style.cssText = "")

}


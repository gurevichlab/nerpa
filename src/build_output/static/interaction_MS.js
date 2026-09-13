// -- spectrum --
let spectrum;
let spectrumNotMatched;
let spectrumCanvas;
// a lookup dictionary: mass value -> experimental_peak_idx
const massPeakIdx = new Map(); 

window.addEventListener("resize", (event) => {
    resizeSpectrum();
});    

function drawSpectrum(spectrumID, variantID){
    document.getElementById('spectrumCanvas').innerHTML = "";
    const svgSpectrum = document.getElementById('spectrumSvg');

    const rect = document.getElementById('spectrum').getBoundingClientRect()
    const canvasWidth = Math.round(rect.width * 0.9);
    const canvasHeight = Math.round(rect.height * 0.7);
    spectrumCanvas = new ChemDoodle.PerspectiveCanvas('spectrumCanvas', canvasWidth, canvasHeight);
    spectrumCanvas.styles.plots_color="grey";
    spectrumCanvas.styles.plots_width= 1;
    spectrumCanvas.styles.text_font_size = 14;
    spectrumCanvas.styles.text_font_families[0] = "Arial";
    spectrumCanvas.styles.text_font_families[1] = "Charcoal";
    spectrumCanvas.styles.text_font_families[2] = 'sans-serif';
    

    const match = getSpectrumVariantMatch(spectrumID, variantID);

    const matchedPeaks = findMatchedPeaks(match, spectrumID);
    const matchedSpectrumJcampFile = peakArrayToJcamp(matchedPeaks);
    spectrum = ChemDoodle.readJCAMP(matchedSpectrumJcampFile); 
   
    spectrumCanvas.loadSpectrum(spectrum);

    const notMatchedPeaks = findNotMatchedPeaks(match, spectrumID);
    const notMatchedSpectrumJcampFile = peakArrayToJcamp(notMatchedPeaks);
    spectrumNotMatched = ChemDoodle.readJCAMP(notMatchedSpectrumJcampFile); 

    spectrumCanvas.loadSpectrum(spectrumNotMatched); 
    colorMatchedPeaks(spectrumCanvas, spectrumID, variantID);

    const oldRepaint = spectrumCanvas.repaint;
    spectrumCanvas.repaint = function(e) {
        oldRepaint.call(this,e);
        updatePeaks(spectrumCanvas);
    };  


    
    svgSpectrum.childNodes.forEach(peakLine => {
        peakLine.addEventListener('click', e => {
            const nid = peakLine.getAttribute('nid').split(',');
            peakLine.setAttribute("clicked", true);
            selectMS(nid);
        });
    });
 
    
    const infoDiv = document.getElementById("spectrumInfo");
    if(spectra[spectrumID].charge != null){
        infoDiv.innerHTML = `<span style="color: #000000; line-height: 2;"><strong>Spectrum ID: ${spectrumID}</strong></span>
                            <br> <span><strong>Precursor:</strong> ${spectra[spectrumID].precursor_mass.toFixed(3)} Da</span>
                             <span style="margin: 0 10px; color: #bbb;">|</span> 
                              <span><strong>MS Score:</strong> ${match.score}</span> 
                               <span style="margin: 0 10px; color: #bbb;">|</span> 
                             <span><strong>Charge:</strong> +${spectra[spectrumID].charge}</span> 
                             `;
    } else {
        infoDiv.innerHTML = `<span style="color: #000000; line-height: 2;"><strong>Spectrum ID: ${spectrumID}</strong></span>
                            <br> <span><strong>Precursor:</strong> ${spectra[spectrumID].precursor_mass.toFixed(3)} Da</span>
                            <span style="margin: 0 10px; color: #bbb;">|</span> 
                             <span><strong>MS Score:</strong> ${match.score}</span> 
                            `;
    }

    initPeakTable(svgSpectrum);
}
function resizeSpectrum(){
    if(spectrumCanvas === undefined) return;
    const rect = document.getElementById('spectrum').getBoundingClientRect()
    const canvasWidth = Math.round(rect.width * 0.9);
    const canvasHeight = Math.round(rect.height * 0.7);
    spectrumCanvas.resize(canvasWidth,canvasHeight);
    spectrumCanvas.loadSpectrum(spectrumNotMatched); 
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
            peak.setAttribute("clicked", true);
            selectMS(nid);
        });
        const peakMzCell = line.insertCell();
        peakMzCell.innerHTML = peak.getAttribute('mass');  

        const peakIntensityCell = line.insertCell();
        peakIntensityCell.innerHTML = peak.getAttribute('intensity');
        
        const chargeCell = line.insertCell();
        chargeCell.innerHTML = peak.getAttribute('charge');  
  
        const peakMassErrorAbsCell = line.insertCell();
        peakMassErrorAbsCell.innerHTML = `${peak.getAttribute('massErrorAbs')}`;;  

        const peakMassErrorRelCell = line.insertCell();
        peakMassErrorRelCell.innerHTML = `${peak.getAttribute('massErrorRel')}`;
    }
}
function peakArrayToJcamp(peaks){
    let jcamp = `##TITLE=${''}\n##DATA TYPE= MASS SPECTRUM\n##XUNITS= m/z\n##YUNITS= relative abundance\n##PEAK TABLE= (XY..XY)`;
    for (const p of peaks) {
        const x = p.mass;
        const y = p.intensity;
        jcamp = jcamp + `${x}, ${y}\n`
    }
    return jcamp + '##END=\n';
}
function findMatchedPeaks(match, spectrumID){
    const resultsPeaks = [];
    for(const peak of match.matched_peaks){
        const matchedPeak = spectra[spectrumID].peaks.find(p => p.mass === peak.experimental_mz)
        if (matchedPeak){
            resultsPeaks.push(matchedPeak);
        }
    }
    return resultsPeaks
}
function findNotMatchedPeaks(match, spectrumID){
    const resultsPeaks = [];
    for(const peak of spectra[spectrumID].peaks){
        const matchedPeak = match.matched_peaks.find(p => p.experimental_mz === peak.mass)
        if (matchedPeak === undefined ){
            resultsPeaks.push(peak);
        }
    }
    return resultsPeaks
}
function colorMatchedPeaks(spectrumCan, spectrumID, variantID){
    const svgSpectrum = document.getElementById('spectrumSvg');
    svgSpectrum.innerHTML = "";

    const variantObjet = getVariantObject(variantID);
    const old_to_new_mon_map = variantObjet.old_to_new_mon_map
    const match = spectra_matching_results.find((entry) => entry.spectrum_id === spectrumID && entry.structure_id === variantID);

    
    const specWidth = spectrumNotMatched.memory.width;
    const specHeight = spectrumNotMatched.memory.height;
    const specOffsetLeft = spectrumNotMatched.memory.offsetLeft;
    const specOffsetBottom = spectrumNotMatched.memory.offsetBottom;
    const specOffsetTop = spectrumNotMatched.memory.offsetTop;

    const origY = spectrumNotMatched.getTransformedY(0,  spectrumCan.styles, specHeight, specOffsetBottom, specOffsetTop);

    for(const peak of spectrum.data){
       
        const coordx = spectrumNotMatched.getTransformedX(peak.x, spectrumCan.styles, specWidth,specOffsetLeft);
        const coordy = spectrumNotMatched.getTransformedY(peak.y,  spectrumCan.styles, specHeight, specOffsetBottom, specOffsetTop);

        const pid = generateID(spectrumID, peak.x);

        const matched_peak = match.matched_peaks.find(peak => peak.experimental_peak_idx === pid); 
        const nid = translateMask(matched_peak.theoretical_fragment_mask, old_to_new_mon_map);
        const mass = parseFloat(peak.x).toFixed(3);
        const intensity = parseFloat(spectra[spectrumID].peaks[matched_peak.experimental_peak_idx].intensity).toFixed(1);
        const massErrorAbs = parseFloat(matched_peak.theoretical_mz - matched_peak.experimental_mz).toFixed(3);
        const massErrorRel = parseFloat((massErrorAbs / matched_peak.experimental_mz ) * 1000000).toFixed(1);
        const charge = matched_peak.charge;
 

        const peakLine = document.createElementNS(svgns, 'line');
        peakLine.setAttribute("id", `${peak.x}`);
        peakLine.setAttribute("nid", `${nid}`);
        peakLine.setAttribute("clicked", false);
        peakLine.setAttribute("mass", `${mass}`);
        peakLine.setAttribute("intensity", `${intensity}`);
        peakLine.setAttribute("massErrorAbs", `${massErrorAbs}`);
        peakLine.setAttribute("massErrorRel", `${massErrorRel}`);
        peakLine.setAttribute("charge", `${charge}`);
        peakLine.setAttribute('x1', coordx);
        peakLine.setAttribute('y1', origY);
        peakLine.setAttribute('x2', coordx);
        peakLine.setAttribute('y2', coordy);
        peakLine.setAttribute('stroke', '#2a6881');
        peakLine.setAttribute('stroke-width', '4');
        svgSpectrum.appendChild(peakLine);
    }
}
function updatePeaks(spectrumCan){
    
    const specWidth = spectrumNotMatched.memory.width;
    const specHeight = spectrumNotMatched.memory.height;
    const specOffsetLeft = spectrumNotMatched.memory.offsetLeft;
    const specOffsetBottom = spectrumNotMatched.memory.offsetBottom;
    const specOffsetTop = spectrumNotMatched.memory.offsetTop;

    const origY = spectrumNotMatched.getTransformedY(0,  spectrumCan.styles, specHeight,specOffsetBottom, specOffsetTop);
    const origX = spectrumNotMatched.getTransformedX(0,  spectrumCan.styles, specWidth, specOffsetLeft);


    for(const peak of spectrum.data){
       
        const coordx = spectrumNotMatched.getTransformedX(peak.x, spectrumCan.styles, specWidth, specOffsetLeft);
        const coordy = spectrumNotMatched.getTransformedY(peak.y,  spectrumCan.styles, specHeight, specOffsetBottom, specOffsetTop);

        if (coordx < origX) continue;
        const peakLine = document.getElementById(`${peak.x}`);

        peakLine.setAttribute('x1', coordx);
        peakLine.setAttribute('y1', origY);
        peakLine.setAttribute('x2', coordx);
        peakLine.setAttribute('y2', coordy);
    }
}
function translateMask(mask, old_to_new_mon_map){
    const deletedIDs = old_to_new_mon_map.filter(t => t[0] != null && t[1] === null).map(t => t[0]);
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
function generateID(spectrumID, peakMass){
    const spectrum = spectra[spectrumID];
    for(const [peakIdx, peakValues] of spectrum.peaks.entries()){
        if(peakValues.mass === peakMass){
            return peakIdx
        }
    }
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

function drawModGraph(nrpID, variant) {
    const variantObj = getVariantObject(variant);

    const graphData_variant = variantObj.new_record;
    const mod_map = variantObj.old_to_new_mon_map;

    const nodesData_mod = [];
    const edgesData_mod = [];
    Object.entries(graphData_variant.monomers).forEach(mon => {
        const coords = getCoords(mon[0], nrpID);
        nodesData_mod.push({
            "id" : mon[0],
            "label" : `${mon[1].name}_${mon[0]}`,
            "color" : getNodeColorNew(mon[0], mod_map),
            "font" : "26",
            "borderWidthSelected" : 4,
            "x": coords.x,
            "y": coords.y
        });
    });
    graphData_variant.monomer_bonds.forEach(bond => {
        edgesData_mod.push({
            "from" :bond[0][0].toString(),
            "to" : bond[0][1].toString(),
            "color" : bond[1][0].atomic_edge.bond_type === 'AMINO' ? 'blue' : 'red',
            "width" : "2",
            "selectionWidth" : "2",
            "hoverWidth" : "2",
            "arrows" : bond[1][0].atomic_edge.bond_type === 'AMINO' ? 'to' : '',
        });
    });  

    const vis_network = buildVisGraph(nodesData_mod, edgesData_mod, "graphModNew", nrpID);
 
    vis_network.addEventListener('click',  e => {
        deselect();
        vis_network.fit();
        vis_network.selectNodes([]);  
    });
    
    const p = document.getElementById("graphModP");
    const variantIDs = extractIDs(variant);
    p.innerHTML =`Modified Monomer Graph <span style="margin: 0 10px; color: #bbb;">|</span> 
                    <span>rank: ${variantIDs.rank}</span>
                    <span style="margin: 0 10px; color: #bbb;">|</span> 
                    <span>number of modifications:${variantIDs.numMods}</span>
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

function getCoords(id, nrpID){
    const node = monomer_graph[nrpID].nodes.find(n => n.id === id)
    return {
        x: node ? node.x : null,
        y: node ? node.y : null
    }
}
//molecule
function getMoleculeVariantData(variantID, nrpID){
    const variantObjet = getVariantObject(variantID);
    const molData = variantObjet.new_record;
    const data = molecule_image[nrpID];
    const aArray = []
    Object.entries(molData.atoms).forEach(a => {
        const coords = getAtomCoords( a[0], data);
        aArray.push({
            "i" : a[0],
            "l" : a[1].name,
            "x": coords.x,
            "y": coords.y,
            "z": 0.0
        });
    });
    const bArray = [];
    Object.entries(molData.atomic_bonds).forEach(([idx, b]) => {
        bArray.push({
            "i" : parseInt(idx),
            "b" : aArray.findIndex(a => parseInt(a.i) === b[0][0]),
            "e" : aArray.findIndex(a => parseInt(a.i) === b[0][1]),
            "o":  parseInt(b[1].arity)
        });
    });
    const monomerDict = {};
    Object.entries(molData.monomers).forEach(m => 
        monomerDict[`${m[1].name}_${m[0]}`] = m[1].atoms
    );
    
    // !! change highlightAtomColors, highlightBonds when colors are included !!

    return {
        a: aArray,
        b: bArray,
        monomers: monomerDict,
        highlightAtomColors: data.highlightAtomColors,
        highlightBonds: data.highlightBonds
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
function getAtomCoords(id, data){
    const atom = data.a.find(a => a.i === id)
    return {
        x: atom ? atom.x : null,
        y: atom ? atom.y : null
    }

}
// -- variant Graph ---
let variantNetwork;
function drawVariantGraph(nrpID, variantID, spectrumID){
    const variantObj = getVariantObject(variantID);

    const graphData_variant = variantObj.new_record;

    const nodesData_variant = [];
    const edgesData_variant = [];
    Object.entries(graphData_variant.monomers).forEach(mon => {
        const coords = getCoords(mon[0], nrpID);
        nodesData_variant.push({
            "id" : mon[0],
            "label" : `${mon[1].name}_${mon[0]}`,
            "color" : getNodeColorVariant(mon[0], nrpID),
            "font" : "26",
            "borderWidthSelected" : 4,
            "x": coords.x,
            "y": coords.y
        });
    });
    graphData_variant.monomer_bonds.forEach(bond => {
        edgesData_variant.push({
            "from" :bond[0][0].toString(),
            "to" : bond[0][1].toString(),
            "color" : bond[1][0].atomic_edge.bond_type === 'AMINO' ? 'blue' : 'red',
            "width" : "2",
            "selectionWidth" : "2",
            "hoverWidth" : "2",
            "arrows" : bond[1][0].atomic_edge.bond_type === 'AMINO' ? 'to' : '',
        });
    });  

    const vis_network = buildVisGraph(nodesData_variant, edgesData_variant, "variantGraphImage", nrpID);
    variantNetwork = vis_network;
    
    variantNetwork.addEventListener('click',  e => {
        if(e.nodes.length === 0){ 
            deselectMS();
            variantNetwork.fit();
        } else {
            selectMS(e.nodes);
        }
    });

    nodeIdLabel.clear();
    for(const entry of nodesData_variant){
        nodeIdLabel.set(entry.id, entry.label);
    }

    displayVarquestMod(spectrumID, variantID);
}
function getNodeColorVariant(id, nrpID){
    const origNode = monomer_graph[nrpID].nodes.find(n => n.id === id);
    if(origNode){
        return origNode.color
    } else {
        return "#cbcbcb"
    }
}
function displayVarquestMod(spectrumID, variantID){
    const match = getSpectrumVariantMatch(spectrumID, variantID);

    let modSum = 0;
    let subText = '';
    for(const mod of match.modifications){
        if(subText != '') subText += ' + ';

        const initMass = mod.initial_monomer_mass;
        const mass_diff = mod.mass_difference;
        modSum += (initMass + mass_diff);

        const modified_node = variantNetwork.nodes.get(mod.monomer_idx.toString());

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

function getSpectrumVariantMatch(spectrumID, variantID){
    return spectra_matching_results.find(
            (entry) => 
                entry.spectrum_id === spectrumID &&
                entry.structure_id === variantID
            );

}

// -- graph/mol --
function switchGraphMol(){
    const graph = document.getElementById('variantGraphImage');
    const molecule = document.getElementById('moleculeImage');
    const showLabelBtn = document.getElementById('showLabelBtn');

    const toMolecule = graph.style.display === 'block';
  
    molecule.style.display = toMolecule ? 'block' : 'none';
    showLabelBtn.style.display = toMolecule ? 'block' : 'none';
    graph.style.display = toMolecule ? 'none' : 'block';

    document.getElementById('switchBtn').textContent = toMolecule ? "Switch to graph view" : "Switch to molecule view";
}

function selectMS(nid){
    let nidString = nid.map(String);
    let nidInt = nid.map(Number);
    const hasInvalidIdx = nidInt.some(Number.isNaN);

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

    // select peaks
    const svgSpectrum = document.getElementById('spectrumSvg');
    const spectrumClick =  Array.from(svgSpectrum.childNodes).some(peakLine => peakLine.getAttribute('clicked') === 'true');
    for(const peak of svgSpectrum.childNodes){
        if(!peak.id) continue;
        const nidPeak = peak.getAttribute('nid').split(',');
        const sameNids = nidPeak.length === nidString.length && nidPeak.every(val => nidString.includes(val));
        const containNids = nidPeak.some(val => nidString.includes(val));
        if(sameNids){
            peak.setAttribute("opacity", "1");
            peak.removeAttribute('stroke-dasharray');
        } else if (containNids && !spectrumClick){
            peak.setAttribute("opacity", "1");
            peak.setAttribute('stroke-dasharray', '2 2');
        } else {
            peak.setAttribute("opacity", "0.3");
            peak.removeAttribute('stroke-dasharray');
        }
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

}
function deselectMS(){

    // deselect variant graph
    variantNetwork.selectNodes([]); 
    increaseTranparency([], variantNetwork.nodes, variantNetwork.edges);

    const svgSpectrum = document.getElementById('spectrumSvg');
    for(const peak of svgSpectrum.childNodes){
        if(!peak.id) continue;
        peak.setAttribute("opacity", "1");
        peak.removeAttribute('stroke-dasharray');
        peak.removeAttribute('clicked');
    }

    // deselect peak Table
    const peakTable = document.querySelector('#peakTable tbody');
    [...peakTable['rows']].forEach(r => r.style.cssText = "")

}


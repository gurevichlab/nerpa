// -- spectrum --
let spectrum;
let spectrumNotMatched;
// a lookup dictionary: mass value -> experimental_peak_idx
const massPeakIdx = new Map(); 

function loadSpectrum(spectrumID, structureID){
    document.getElementById('spectrumCanvas').innerHTML = "";
    const svgSpectrum = document.getElementById('spectrumSvg');

    const rect = svgSpectrum.getBoundingClientRect();
    let spectrumCan = new ChemDoodle.PerspectiveCanvas('spectrumCanvas', 600, 400);
    spectrumCan.styles.plots_color="grey";
    spectrumCan.styles.plots_width= 1;
    spectrumCan.styles.text_font_size = 14;
    spectrumCan.styles.text_font_families[0] = "Arial";
    spectrumCan.styles.text_font_families[1] = "Charcoal";
    spectrumCan.styles.text_font_families[2] = 'sans-serif';
    

    const matchedPeaks = findMatchedPeaks(spectrumID, structureID);
    const matchedSpectrumJcampFile = peakArrayToJcamp(matchedPeaks)
    spectrum = ChemDoodle.readJCAMP(matchedSpectrumJcampFile); 
   
    spectrumCan.loadSpectrum(spectrum);

    const notMatchedPeaks = findNotMatchedPeaks(spectrumID, structureID);
    const notMatchedSpectrumJcampFile = peakArrayToJcamp(notMatchedPeaks)
    spectrumNotMatched = ChemDoodle.readJCAMP(notMatchedSpectrumJcampFile); 

    spectrumCan.loadSpectrum(spectrumNotMatched); 
    highlightPeaks(spectrumCan, spectrumID, structureID);

    const oldRepaint = spectrumCan.repaint;
    spectrumCan.repaint = function(e) {
        oldRepaint.call(this,e);
        updatePeaks(spectrumCan);
    };  
    
    svgSpectrum.childNodes.forEach(peakLine => {
        peakLine.addEventListener('click', e => {
            const nid = peakLine.getAttribute('nid').split(',');
            peakLine.setAttribute("clicked", true);
            select(nid,extractNrpID(structureID));
        });
    });
 
    
    const tSpectra = document.querySelector("#spectrumInfo table");
    tSpectra.innerHTML = '';
    const tHead = tSpectra.createTHead();
    const hLine = document.createElement("tr");
    const thName = document.createElement("th");
    thName.textContent = 'Spectrum name';
    hLine.appendChild(thName);
    const thMass = document.createElement("th");
    thMass.textContent = 'Precursor mass';
    hLine.appendChild(thMass);
    if(spectra[spectrumID].charge != null){
        const thCharge = document.createElement("th");
        thCharge.textContent = 'Charge';
        hLine.appendChild(thCharge);
    }
    tHead.appendChild(hLine);
    const line = tSpectra.insertRow();
    if(spectra[spectrumID].charge != null){
        const spCharge = line.insertCell(0);
        spCharge.innerHTML = spectra[spectrumID].charge; 
    }
    const spPreMass = line.insertCell(0);
    spPreMass.innerHTML = spectra[spectrumID].precursor_mass.toFixed(3);
    const spName = line.insertCell(0);
    spName.innerHTML = spectrumID;
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

function findMatchedPeaks(spectrumID, structureID){
    const match = spectra_matching_results.find((entry) => entry.spectrum_id === spectrumID && entry.structure_id === structureID);
    const resultsPeaks = [];
    for(const peak of match.matched_peaks){
        const matchedPeak = spectra[spectrumID].peaks.find(p => p.mass === peak.experimental_mz)
        if (matchedPeak ){
            resultsPeaks.push(matchedPeak);
        }
    }
    return resultsPeaks
}

function findNotMatchedPeaks(spectrumID, structureID){
    const match = spectra_matching_results.find((entry) => entry.spectrum_id === spectrumID && entry.structure_id === structureID);
    const resultsPeaks = [];
    for(const peak of spectra[spectrumID].peaks){
        const matchedPeak = match.matched_peaks.find(p => p.experimental_mz === peak.mass)
        if (matchedPeak === undefined ){
            resultsPeaks.push(peak);
        }
    }
    return resultsPeaks
}

function highlightPeaks(spectrumCan, spectrumID, structureID){
    const svgSpectrum = document.getElementById('spectrumSvg');
    svgSpectrum.innerHTML = "";
    
    const specWidth = spectrumNotMatched.memory.width;
    const specHeight = spectrumNotMatched.memory.height;
    const specOffsetLeft = spectrumNotMatched.memory.offsetLeft;
    const specOffsetBottom = spectrumNotMatched.memory.offsetBottom;
    const specOffsetTop = spectrumNotMatched.memory.offsetTop;

    const origX = spectrumNotMatched.getTransformedX(0, spectrumCan.styles, specWidth,specOffsetLeft);
    const origY = spectrumNotMatched.getTransformedY(0,  spectrumCan.styles, specHeight,specOffsetBottom, specOffsetTop);

    for(const peak of spectrum.data){
       
        const coordx = spectrumNotMatched.getTransformedX(peak.x, spectrumCan.styles, specWidth,specOffsetLeft);
        const coordy = spectrumNotMatched.getTransformedY(peak.y,  spectrumCan.styles, specHeight, specOffsetBottom, specOffsetTop);

        const pid = generateID(spectrumID, peak.x);

        const match = spectra_matching_results.find((entry) => entry.spectrum_id === spectrumID && entry.structure_id === structureID);
        const matched_peak = match.matched_peaks.find(peak => peak.experimental_peak_idx === pid); 
        const nid = translateMask(matched_peak.theoretical_fragment_mask);
        const mass = parseFloat(peak.x).toFixed(3);
        const intensity = parseFloat(spectra[spectrumID].peaks[matched_peak.experimental_peak_idx].intensity).toFixed(1);
        const massErrorAbs = parseFloat(matched_peak.theoretical_mz - matched_peak.experimental_mz).toFixed(3);
        const massErrorRel = parseFloat((massErrorAbs / matched_peak.experimental_mz ) * 1000000).toFixed(1);

        const peakLine = document.createElementNS(svgns, 'line');
        peakLine.setAttribute("id", `${peak.x}`);
        peakLine.setAttribute("nid", `${nid}`);
        peakLine.setAttribute("clicked", false);
        peakLine.setAttribute("mass", `${mass}`);
        peakLine.setAttribute("intensity", `${intensity}`);
        peakLine.setAttribute("massErrorAbs", `${massErrorAbs}`);
        peakLine.setAttribute("massErrorRel", `${massErrorRel}`);
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

function translateMask(mask){
    const results = [];
    let id = 1;
    for(const char of mask.split('')){
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
function drawModGraph(nrpID) {
    const graphMod = Object.values(candidate_NRPs).find(entry => entry.compound_id === nrpID);
    const graphData_new = Object.values(graphMod.new_variants)[1].new_record;
    const nodesData_new = [];
    const edgesData_new = [];
    Object.entries(graphData_new.monomers).forEach(mon => {
        const coords = getCoords(mon[0]);
        nodesData_new.push({
            "id" : mon[0],
            "label" : `${mon[1].name}_${mon[0]}`,
            "color" : "#cbcbcb",
            "font" : "26",
            "borderWidthSelected" : 4,
            "x": coords.x,
            "y": coords.y
        });
    });
    graphData_new.monomer_bonds.forEach(bond => {
        edgesData_new.push({
            "from" :bond[0][0].toString(),
            "to" : bond[0][1].toString(),
            "color" : "blue",
            "width" : "2",
            "selectionWidth" : "2",
            "hoverWidth" : "2",
            "arrows" : "2",
        });
    });  

    const graphData_old = graphMod.original_record;
    const nodesData_old = [];
    const edgesData_old = [];
    Object.entries(graphData_old.monomers).forEach(mon => {
        const coords = getCoords(mon[0]);
        nodesData_old.push({
            "id" : mon[0],
            "label" : `${mon[1].name}_${mon[0]}`,
            "color" : "#cbcbcb",
            "font" : "26",
            "borderWidthSelected" : 4,
            "x": coords.x,
            "y": coords.y
        });
    });
    graphData_old.monomer_bonds.forEach(bond => {
        edgesData_old.push({
            "from" :bond[0][0].toString(),
            "to" : bond[0][1].toString(),
            "color" : "blue",
            "width" : "2",
            "selectionWidth" : "2",
            "hoverWidth" : "2",
            "arrows" : "2",
        });
    });  


    const nodes_old = new vis.DataSet(
                nodesData_old.map(n => ({...n}))
            );
    const edges_old = new vis.DataSet(
                edgesData_old.map(e => ({...e}))
            );
    const nodes_new = new vis.DataSet(
                nodesData_new.map(n => ({...n}))
            );
    const edges_new = new vis.DataSet(
                edgesData_new.map(e => ({...e}))
            );
    const graph_div_old = document.getElementById("graphModOld");
    const graph_div_new = document.getElementById("graphModNew");
    const graph_data_old = { nodes: nodes_old, edges: edges_old };
    const graph_data_new = { nodes: nodes_new, edges: edges_new };
    const graph_options = { 
        physics: { enabled: false},
        edges: {
            arrows: {
                to: {enabled: true, scaleFactor: 1.5, type: 'arrow'},
            },
            smooth: {
                enabled: true,
                type: 'dynamic', 
                roundness: 0.5      
            }
        },
        interaction:{
            zoomView:false,
            multiselect: true,
        },
    };

    const networkOld = new vis.Network(graph_div_old, graph_data_old, graph_options);
    const networkNew = new vis.Network(graph_div_new, graph_data_new, graph_options);
    Object.values(graphMod.new_variants)[1].old_to_new_mon_map.forEach(tup => {
        if(tup[1] === null){
            //(old_idx, null) -- the monomer old_idx was removed
            nodes_old.update({
                id: tup[0].toString(),
                color: "#e79898"    
            });
        } else if(tup[0] === null){
            //(null, new_idx) -- the monomer new_idx was inserted
            nodes_new.update({
                id: tup[1].toString(),
                color: "#98e799" 
            });
        } else {
            //(old_idx, new_idx) -- the monomer old_idx was substituted (in this case old_idx=new_idx)
            nodes_old.update({
                id: tup[0].toString(),
                color: "#ffff9e"    
            });
            nodes_new.update({
                id: tup[1].toString(),
                color: "#ffff9e"    
            });

        }
    });

}

function getCoords(id){
    const node = nodes.get(id);
    return {
        x: node ? node.x : null,
        y: node ? node.y : null
    }
}

// -- varQuest modification --
function displayVarquestMod(spectrumID, structureID){
    const match = spectra_matching_results.find((entry) => entry.spectrum_id === spectrumID && entry.structure_id === structureID);

    let subText = '';
    for(const mod of match.modifications){
        if(subText != '') subText += ' + ';
        const initMass = mod.initial_monomer_mass.toFixed(3);
        const mass_diff = mod.mass_difference.toFixed(3);
        const modMonNode = nodes.get(mod.monomer_idx.toString())
        const color = modMonNode.color;
        const border = darkenColor(color, 30);
        if (mass_diff < 0){
            subText += `<span style="
                            background-color: ${color};
                            border: 2px dashed ${border};
                            border-radius: 999px;
                            padding: 2px 6px;
                            ">${modMonNode.label}</span>
                             - ${mass_diff * -1} Da`;
        } else {
            subText += `<span style="background-color: ${color};
                            border: 2px dashed ${border};
                            border-radius: 999px;
                            padding: 2px 6px;">${modMonNode.label}</span>
                                     + ${mass_diff} Da`;
        }
        nodes.update({
            id: modMonNode.id,
            shapeProperties: {
                borderDashes: [5, 5] 
            },
            "borderWidth": 2,    
        });
        nodes.add({
            id: Math.PI,
            label: `${mass_diff} Da`,
            x: modMonNode.x + 60,
            y: modMonNode.y + 30,
            "fixed": true,
            "physics": false, 
            "labelHighlightBold": false,
            "color": {
                "background":  "transparent", 
                "border":  "transparent",    
                "highlight": {
                    "background":  "transparent", 
                    "border":  "transparent",   
                }
            },
                                        
            font: { 
                color: '#444',
                size: 16, 
                mod: 'bold'
            }
        }); 
    }

    const spectrumMass = match.spectrum_mass.toFixed(3);

    const textHtml = `${subText} = ${spectrumMass} Da`;
   

    const gdiv = document.getElementById('modTable');
    gdiv.innerHTML = '';

    const tMod = document.createElement("table");
  

    const tHead = tMod.createTHead();
    const hLine = document.createElement("tr");
    const thLabel = document.createElement("th");
    thLabel.textContent = 'MS-based mod:';
    hLine.appendChild(thLabel);
    tHead.appendChild(hLine);
    const line = tMod.insertRow();
    const lineContent = line.insertCell(0);
    lineContent.innerHTML = textHtml;  
    gdiv.appendChild(tMod);


}

function darkenColor(hex, percent) {
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
// -- graph/mol --

function switchGraphMol(){
    const graph = document.getElementById('graphImage');
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
            const peakDiv = document.getElementById("peak");
            peakDiv.style.display = 'block';
            document.getElementById("peakMz").textContent = peak.getAttribute('mass');
            document.getElementById("peakIntensity").textContent = peak.getAttribute('intensity');
            document.getElementById("peakMassErrorAbs").textContent = `${peak.getAttribute('massErrorAbs')} Da`;
            document.getElementById("peakMassErrorRel").textContent = `${peak.getAttribute('massErrorRel')} ppm`;
        } else if (containNids && !spectrumClick){
            peak.setAttribute("opacity", "1");
            peak.setAttribute('stroke-dasharray', '2 2');
        } else {
            peak.setAttribute("opacity", "0.3");
            peak.removeAttribute('stroke-dasharray');
        }
        peak.setAttribute('clicked', false);
    }

}
function deselectMS(){
    const peakDiv = document.getElementById("peak");
    const svgSpectrum = document.getElementById('spectrumSvg');
    peakDiv.style.display = 'none';
    for(const peak of svgSpectrum.childNodes){
        if(!peak.id) continue;
        peak.setAttribute("opacity", "1");
        peak.removeAttribute('stroke-dasharray');
        peak.removeAttribute('clicked');
    }

}


// -- local variables --
const svgns = "http://www.w3.org/2000/svg";
// --- Module Viewer ---
function drawModuleViewerAlign(bgcID, nrpID){
    const svg = document.getElementById("moduleViewerAlign");

    const bgcKey = bgcID.constructor == Object ? makeBgcKey(bgcID) : bgcID;
    const match = data.find(item =>
                            makeBgcKey(item["bgc_variant_id"]["bgc_id"]) === bgcKey &&
                            item["nrp_variant_id"]["nrp_id"] === nrpID);
    
    // this is necessary because match.alignments is currently saved as a list within a list of length 1
    const array_in_array =  Array.isArray(match.alignments) && 
                            match.alignments.length === 1 &&
                            Array.isArray(match.alignments[0])
    const alignments = array_in_array ? match.alignments[0] : match.alignments; 
         
    const genes = prepareGenesAlign(bgcKey, alignments);

    const drawData = buildDrawData(genes);

    const lastGene = drawData[drawData.length - 1];
    const viewerWidth = lastGene.position.x + lastGene.length + 15;
    svg.setAttribute("width", viewerWidth);
    svg.innerHTML = "";

    drawModuleViewer(drawData, svg);
    addClickHandler(svg);
}

function drawModuleViewerBGC(bgcID, nrpID){
    const svgBGC = document.getElementById("moduleViewerBGC");
    
    const bgcKey = bgcID.constructor == Object ? makeBgcKey(bgcID) : bgcID;
    const match = data.find(item =>
                            makeBgcKey(item["bgc_variant_id"]["bgc_id"]) === bgcKey &&
                            item["nrp_variant_id"]["nrp_id"] === nrpID);
     
    // this is necessary because match.alignments is currently saved as a list within a list of lenth 1
    const array_in_array =  Array.isArray(match.alignments) && 
                            match.alignments.length === 1 &&
                            Array.isArray(match.alignments[0])
    const alignments = array_in_array ? match.alignments[0] : match.alignments;
    const genes = prepareGenesBGC(bgcKey, alignments);

    const drawData = buildDrawData(genes, true);

    const lastGene = drawData[drawData.length - 1];
    const viewerWidth = lastGene.position.x + lastGene.length + 15;
    svgBGC.setAttribute("width", viewerWidth);
    svgBGC.innerHTML = "";

    drawModuleViewer(drawData, svgBGC);
    addClickHandler(svgBGC);
}
function prepareGenesAlign(bgcKey, alignments){
    const genesMatch = structuredClone(bgcModuleIndex.get(bgcKey));
    const genes = [];
    let pendingAlignments = undefined;
    let lastGeneID = undefined;

    // assign each module its corresponding alignment 
    for(const al of alignments){
        // no A domain idx -> no module
        if(al["A-domain_idx"] === "---") {
            continue;
        }
        // init pendingAlignmets and lastGeneID
        if(pendingAlignments === undefined) {
            pendingAlignments = [al];
            lastGeneID = al.Gene;
            continue;
        }
        
        if(al["A-domain_idx"] === 0){
            const g = structuredClone(genesMatch.find(item => item.gene_id === lastGeneID));
            let idx = 0;
            for(const m of g.modules){
                if(m.a_domain === null) {
                    m.alignment = undefined;
                } else {
                    m.alignment = pendingAlignments[idx];
                    idx ++;
                }
            }
            genes.push(g);

            // reinit pendingAlignmets and lastGeneID
            pendingAlignments = [al];
            lastGeneID = al.Gene;

        } else {
            pendingAlignments.push(al);
        }
    }

    // assign alignmnet to modules of last gene 
    const g = genesMatch.find(item => item.gene_id === lastGeneID);
    let idx = 0;
    for(const m of g.modules){
        if(m.a_domain === null) {
            m.alignment = undefined;
        } else {
            m.alignment = pendingAlignments[idx];
            idx ++;
        }
    }
    genes.push(g);

    return genes;
}
function prepareGenesBGC(bgcKey, alignments){
    const genes = structuredClone(bgcModuleIndex.get(bgcKey));

    // assign each module its corresponding alignment
    for(const g of genes){
        let idx = 0;
        for(const m of g.modules){
            if(m.a_domain === null) {
                m.alignment = [];
            } else {
                m.alignment = alignments.filter(
                    item =>
                            item.Gene === g.gene_id &&
                            item["A-domain_idx"] === idx
                );
                idx ++;
            }
        }
    }
    return genes
}
function buildDrawData(geneList, r=false){
    /* gene object
    {
        id: Number,
        postion:  {x: Number, y: Number},
        length: Number (in px),
        displayed: Bool,
        nid: Number, (if of graph nodes matched to  all modules)
        reversed: Bool,
        modules: Array<Module>
    }
    module object
    {
        id: Number,
        nid: Number, (id of graph nodes and molecules)
        domains: Array<String>,
        residue: String,
        displayed: Boolean,
        postion:  {x: Number, y: Number},
        length: Number,
        rad: Number,
    }
    */
    const radFixed = 20; // radius of domain circles
    const dist = 15;   // distance between modules
    const yFixed = 75;

    let drawData = []; // 

    for(const g of geneList){ 
        let gene ={
            id: g.gene_id,                         
            position: {x: undefined, y: yFixed},                          
            length: undefined,  
            displayed: true,  
            nid: undefined,
            reversed: g.coords.strand === 'REVERSE' && r,
            modules: g.modules.map((m, idx)=> {
                return{
                    id: idx +1,
                    nid: m["a_domain"] != null ? 
                                (r ? m.alignment.map(al => al["rBAN_idx"]) : m.alignment["rBAN_idx"]): [''],
                    domains: m.domains_sequence, 
                    residue: selectResidue(m, r),
                    displayed: m["a_domain"] != null,        
                    position: {x: undefined, y: yFixed}, 
                    length: m.domains_sequence.length * radFixed * 2,     
                    rad: radFixed,
                    skipped: r ? false : !Number.isInteger(m.alignment.rBAN_idx),               
                };
            }),              
        };
        gene.length = gene.modules.reduce((c,m) => c + m.length + dist, -dist);
        drawData.push(gene);
    }
    // calc positions of genes and modules
    let xGene = 0;
    for(const gObj of drawData){
        gObj.position.x = xGene + dist;
        xGene = xGene + dist + gObj.length;

        let xModule = gObj.position.x;
        if(gObj.reversed) gObj.modules.reverse();
        let allNotDisplayed = true;
        for(const mObj of gObj.modules){
            if(mObj.displayed) allNotDisplayed = false;
            mObj.position.x = xModule;
            xModule = xModule + dist + mObj.length;
        }
        gObj.displayed = !allNotDisplayed;
        gObj.nid = gObj.modules.map(m => m.nid).flat();
    }

    return drawData;
}
function drawModuleViewer(drawData, svg){
    for(const g of drawData){
        for (const m of g.modules){
            drawModule(g.id, m.id, m.nid, m.position.x, m.position.y, m.length, m.domains, m.residue, m.displayed, m.skipped, m.rad, svg);
        }
        drawGeneArrow(g.id, g.nid, g.position.x, g.position.y, g.length, g.displayed, g.reversed, svg);
    }
}

function drawModule(gid, mid, nid, x, y, length, domains, residue, displayed, skipped, rad, svg){
        
    const module = document.createElementNS(svgns, 'g');
    module.setAttribute("id", `M${mid},G${gid},${nid}`);
      
    let i = 0;
    let xi;
    for(const d of domains){
        xi = rad + i * (rad * 2) + x;
        i++;
        const circle = document.createElementNS(svgns, 'circle');
        let y_circ = y;
        let domainName;
        let typeIdx = '';
        switch(d) {
            case 'A':
                circle.setAttribute('fill', `rgb(${179}, ${130}, ${238})`); 
                domainName = 'A';
                break;
            case 'PCP':
                circle.setAttribute('fill', `rgb(${142}, ${189}, ${242})`); 
                domainName = 'CP';
                break;
            case 'E':
                circle.setAttribute('fill', `rgb(${127}, ${127}, ${238})`);  
                domainName = 'E';
                break;
            case 'PKS':
                circle.setAttribute('fill', `rgb(${110}, ${248}, ${126})`); 
                domainName = 'KS';
                break;
            // C domain
            case 'C_STARTER':
                typeIdx = 'St';
                circle.setAttribute('fill', `rgb(${128}, ${128}, ${240})`); 
                domainName = 'C';
                break;
            case 'C_LCL':
                typeIdx = 'LCL';
                circle.setAttribute('fill', `rgb(${128}, ${128}, ${240})`); 
                domainName = 'C';
                break;
            case 'C_DCL':
                typeIdx = 'DCL';
                circle.setAttribute('fill', `rgb(${128}, ${128}, ${240})`); 
                domainName = 'C';
                break;
            case 'C_DUAL':
                typeIdx = 'DUAL';
                circle.setAttribute('fill', `rgb(${128}, ${128}, ${240})`); 
                domainName = 'C';
                break; 
            case 'C':
                circle.setAttribute('fill', `rgb(${128}, ${128}, ${240})`); 
                domainName = 'C';
                break;
            case 'TE_TD':
                circle.setAttribute('fill', `rgb(${245}, ${187}, ${238})`);
                domainName = 'TE'; 
                break;
            case 'TE':
                circle.setAttribute('fill', `rgb(${245}, ${187}, ${238})`);
                domainName = 'TE'; 
                break;
            case 'TD':
                circle.setAttribute('fill', `rgb(${245}, ${187}, ${238})`);
                domainName = 'TD'; 
                break;
            case 'MT':
                circle.setAttribute('fill', `rgb(${213}, ${213}, ${213})`); 
                y_circ = y - 12;
                domainName = 'nMT';
                break;
            case 'CTERM':
                // color? 
                break;
            case 'NTERM':
                // color? 
                break;
            default:
                throw new Error(`'${d}'' is no domain.`);
        }
        circle.setAttribute('cx', xi);
        circle.setAttribute('cy', y_circ);
        circle.setAttribute('r', rad); 
        circle.setAttribute('stroke', 'black');
        module.appendChild(circle);
        // domain text
        const letter = document.createElementNS(svgns, 'text');
        letter.setAttribute('x', xi);
        letter.setAttribute('y', y_circ );
        letter.setAttribute('font-family', 'Arial');
        letter.setAttribute('font-size', '16');
        letter.setAttribute('fill', 'black');
        letter.setAttribute('text-anchor', 'middle');
        letter.setAttribute('dominant-baseline', 'central');
        letter.textContent = domainName;
        const type = document.createElementNS(svgns, "tspan");
        type.textContent = typeIdx;
        type.setAttribute("baseline-shift", "sub");
        type.setAttribute("font-size", "60%");
        letter.appendChild(type);
        module.appendChild(letter);
    }
    // line below domains
    const line = document.createElementNS(svgns, 'line');
    line.setAttribute('x1', x);
    line.setAttribute('y1', y + 30);
    line.setAttribute('x2', length + x);
    line.setAttribute('y2', y + 30);
    line.setAttribute('stroke', 'black');
    line.setAttribute('stroke-width', '1.5');
    module.appendChild(line);
    // id + residue description
    const idMod = document.createElementNS(svgns, 'text');
    idMod.setAttribute('x', length / 2 + x);
    idMod.setAttribute('y', y + 40);
    idMod.setAttribute('font-family', 'Arial');
    idMod.setAttribute('font-size', '16');
    idMod.setAttribute('fill', 'black');
    idMod.textContent = `M${mid}`;
    idMod.setAttribute('text-anchor', 'middle');
    idMod.setAttribute('dominant-baseline', 'central');
    module.appendChild(idMod);
    const res = document.createElementNS(svgns, 'text');
    res.setAttribute('x', length / 2 + x);
    res.setAttribute('y', y + 56);
    res.setAttribute('font-family', 'Arial');
    res.setAttribute('font-size', '15');
    res.setAttribute('fill',  `rgb(${85}, ${85}, ${85})`);
    res.setAttribute('text-anchor', 'middle');
    res.setAttribute('dominant-baseline', 'central');
    res.textContent = residue;
    module.appendChild(res);

    // rectangles
    if(skipped){
        const bracketL = document.createElementNS(svgns, 'path');
        bracketL.setAttribute("d", `M ${x + 5} ${y + 61} L${x - 3} ${y + 61} L${x - 3} ${y - 30} L${x + 5} ${y - 30}`);
        bracketL.setAttribute('fill', 'none');
        bracketL.setAttribute('stroke', '#cc5253');
        bracketL.setAttribute('stroke-width', '1');
        module.appendChild(bracketL);
        const bracketR = document.createElementNS(svgns, 'path');
        bracketR.setAttribute("d", `M ${x + length - 5} ${y + 61} L${x + length + 3} ${y + 61} L${x + length + 3} ${y - 30} L${x + length - 5} ${y - 30}`);
        bracketR.setAttribute('fill', 'none');
        bracketR.setAttribute('stroke', '#cc5253');
        bracketR.setAttribute('stroke-width', '1');
        module.appendChild(bracketR); 
        const skipWord = document.createElementNS(svgns, 'text');
        skipWord.setAttribute('x', length / 2 + x);
        skipWord.setAttribute('y', y - 30);
        skipWord.setAttribute('font-family', 'Arial');
        skipWord.setAttribute("font-size", "60%");
        skipWord.setAttribute('fill', '#cc5253');
        skipWord.setAttribute('text-anchor', 'middle');
        skipWord.setAttribute('dominant-baseline', 'central');
        skipWord.textContent = 'skipped';
        module.appendChild(skipWord); 


    }
    
    if(!displayed) {
        module.setAttribute("opacity", "0.3");
        module.setAttribute("displayed", "false");
    }
    module.displayed = displayed;
    svg.appendChild(module);
}

function drawGeneArrow(id, nid, x, y, length, displayed, reversed, svg){
    const geneEnd = x + length;
    const geneStart = x;
    let arrowStart = x;
    let arrowhead = x + length;
    let arrowlineEnd = x + (length * 0.9);
    const geneArr = document.createElementNS(svgns, 'g');
    geneArr.setAttribute("id", `G${id},${nid}`);

    // draw arrow
    const arrow = document.createElementNS(svgns, 'path');
    if(reversed){
        arrowStart = geneEnd;
        arrowhead = geneStart;
        arrowlineEnd = geneEnd - (length * 0.9);
    }
    let d = ` M ${arrowStart} ${y - 60} L ${arrowlineEnd} ${y - 60} L ${arrowlineEnd} ${y - 65} L ${arrowhead} ${y - 50} L ${arrowlineEnd} ${y - 35} L ${arrowlineEnd} ${y - 40} L ${arrowStart} ${y - 40} Z`;
    arrow.setAttribute("d", d);
    arrow.setAttribute('fill', 'black');
    geneArr.appendChild(arrow);

    // set arrow description
    const arrowDes = document.createElementNS(svgns, 'text');
    arrowDes.setAttribute('x',  geneStart + (length / 2));
 arrowDes.setAttribute('y', y - 50);
    arrowDes.setAttribute('font-family', 'Arial');
    arrowDes.setAttribute('font-size', '14');
    arrowDes.setAttribute('fill', 'white');
    arrowDes.setAttribute('text-anchor', 'middle');
    arrowDes.setAttribute('dominant-baseline', 'central');
    if(id.length * 7 > length) {
        const separator = "...";
        const available = (length / 7 ) - separator.length;
        const prefixLength = Math.floor(available / 2);
        const suffixLength = available  - prefixLength;
        const preffix = id.slice(0, prefixLength) 
        const suffix =  id.slice(-suffixLength);
        arrowDes.textContent = preffix + separator + suffix;
    }
    else { 
        arrowDes.textContent = id;
    }
    geneArr.appendChild(arrowDes);

    if(!displayed) {
        geneArr.setAttribute("opacity", "0.3");
        geneArr.setAttribute("displayed", "false");
    }
    geneArr.displayed = displayed;
    svg.appendChild(geneArr);
}

function splitId(id){
    const split = id.split(',');

    const sDict = {};
    sDict.nid = [];
    for(const s of split){
        if(s.startsWith("M")) sDict.mid = s;
        else if(s.startsWith("G")) sDict.gid = s;
        else sDict.nid.push(s);
    }
    
    return sDict;
}

function selectResidue(m, bgc){
    // module in BGC order can have multiple alignments but all have same 'Top_scoring_residues' 
    const alignment = Array.isArray(m.alignment) ? m.alignment[0] : m.alignment;
    if(bgc){
        return (
            m["a_domain"] != null ?
                ( 
                    alignment["Top_scoring_residues"] != "---" && alignment["Top_scoring_residues"] != "_UNKNOWN"  ?
                    ( 
                        m["domains_sequence"].includes('E') ?
                            `D-${alignment["Top_scoring_residues"]}`
                        : alignment["Top_scoring_residues"] 
                    )
                    : 'unk'
                )
            : 'x'  
        )
    } else {
        return (
            m["a_domain"] != null ?
                (
                    alignment["Top_scoring_residues"] != "---" && alignment["Top_scoring_residues"] != "_UNKNOWN"  ?
                    (
                        alignment["Modifying_domains"] ===  "---" ?
                            alignment["Top_scoring_residues"] 
                        : `D-${alignment["Top_scoring_residues"]}`
                    )
                    : 'unk'
                )
            : 'x'  
        )
    }
}

function addClickHandler(svg){
    svg.childNodes.forEach(c => {
        const nid = splitId(c.id).nid.filter(n => n != '');
        if(nid.length > 0){
            c.addEventListener('click', e => {
                c.setAttribute('clicked', 'true');
                select(nid);
                c.removeAttribute('clicked');
            });
        }
    });
}
// --- Graph ---

// -- graph variables -- 
let network;
let nodes;
let edges; 
let hasSelectedNodes = false;
// a lookup dictionary: node ID -> node label
const nodeIdLabel = new Map(); 

function drawGraph(nrpID){
    const graph = monomer_graph[nrpID];
    const graph_div = document.getElementById("graphImage");
    const [graph_data, graph_options] = buildGraph(graph);
    network = new vis.Network(graph_div, graph_data, graph_options);
    
    network.addEventListener('click',  e => {
        if(e.nodes.length === 0){ 
            deselect();
            network.fit();
        } else {
            select(e.nodes);
        }
    });

    maxZoomGraph();
    
    const idsToUpdate = chiralityCheck(nrpID);
    for(const id of idsToUpdate){
        const node = nodes.get(id);
        node.label = `D-${node.label}`;
    }
    nodeIdLabel.clear();

    for(const entry of graph.nodes){
        nodeIdLabel.set(entry.id, entry.label);
    }
}

function buildGraph(data) {
    
    nodes = new vis.DataSet(
                data.nodes.map(n => ({...n}))
            );
    edges = new vis.DataSet(
                data.edges.map(e => ({...e}))
            );
    
    const graph_data = { nodes: nodes, edges: edges };
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
            },
            arrowStrikethrough: false,
        },
        interaction:{
            zoomView:true,
            multiselect: true,
        },
    };
 
    return [graph_data, graph_options]

}

function maxZoomGraph(){
    // this code is from https://github.com/visjs/vis-network/issues/574 as vis-network has no minZoom 
    let lastPosition = null;
    const max_zoom = 2;
    const min_zoom = 0.3;
    network.on("zoom", function (params) {
        if (params.scale < min_zoom|| params.scale > max_zoom) { // adjust this value according to your requirement
            network.moveTo({
            position: lastPosition, // use the last position before zoom limit
            scale: params.scale > max_zoom ? max_zoom : min_zoom // this scale prevents zooming out beyond the desired limit
            });
        } else {
            // store the current position as the last position before zoom limit
            lastPosition = network.getViewPosition();
        }
        });
        // on pan, store the current position
        network.on("dragEnd", function () {
        lastPosition = network.getViewPosition();
    });
}
function increaseTranparency(nodeIDs, networkNodes = nodes, networkEdges = edges){
    const realNodes = networkNodes.get({filter: n => n.id >= 0})
    // reset
    realNodes.forEach(n => networkNodes.update({ id: n.id, opacity: 1, font: '26 arial black'}));
    networkEdges.get().forEach(e => networkEdges.update({id: e.id, color: e.color === '#E0A59D' || e.color === 'red' ?'red' : 'blue'}));
    // color graph transparent expect selected node
    if(!(nodeIDs.length === 0)){
        realNodes.forEach(n => !nodeIDs.includes(n.id) && n.id >= 0 ? networkNodes.update({ id: n.id, opacity: 0.3, font: '26 arial #D3D3D3'}) : '');
        networkEdges.get().forEach(e => networkEdges.update({id: e.id, color: e.color === 'red' ? '#E0A59D' : '#94ACD4'}));
    }
}

function chiralityCheck(nrpID){
    let alignments = data.find(item => item.nrp_variant_id.nrp_id === nrpID).alignments[0];

    // this is necessary because match.alignments is currently saved as a list within a list of lenth 1
    const array_in_array =  Array.isArray(alignments) && 
                            alignments.length === 1 &&
                            Array.isArray(alignments[0])
    alignments = array_in_array ? alignments[0] : alignments; 

    const idsToUpdate = alignments.filter(al =>  al.NRP_chirality === 'D').map(item => item.rBAN_idx.toString());
    
    return idsToUpdate
}

// --- Molecule ---
// -- molecule variables --
let molCanvas;

window.addEventListener("resize", (event) => {
        resizeMol();
});    
function drawMolecule(data){
    document.getElementById('moleculeImageCanvas').innerHTML = "";
    const rect = document.getElementById('moleculeImage').getBoundingClientRect();
    molCanvas = new ChemDoodle.ViewerCanvas('moleculeImageCanvas', Math.round(rect.width), Math.round(rect.height));
    
    const mol = new ChemDoodle.io.JSONInterpreter().molFrom(data);
    mol.selected = false;
    mol.highlightBonds = data.highlightBonds;
    mol.highlightAtomColors = data.highlightAtomColors;
    mol.monomers = data.monomers;

    for(const atom of mol.atoms){
        if(atom.label.includes('_')) {
            atom.altLabel = atom.label;
            atom.residue = atom.label;
            atom.label = 'C';
        }    
    }
  
    molCanvas.styles.bonds_width_2D = .9;
    molCanvas.styles.bonds_saturationWidthAbs_2D = 2.6;
    molCanvas.styles.bonds_hashSpacing_2D = 2.5;
    molCanvas.styles.atoms_font_size_2D = 10;
    molCanvas.styles.atoms_font_families_2D = ['Helvetica', 'Arial', 'sans-serif'];
    molCanvas.styles.backgroundColor = 'transparent';

    const bondLen = getMaxBondLen(mol, molCanvas.width, molCanvas.height);
    mol.scaleToAverageBondLength(bondLen);

    molCanvas.loadMolecule(mol);
    paintMolecule(mol);

    enableZoomMol(document.getElementById('moleculeImageCanvas'),mol);    
}

function resizeMol(){
    if(molCanvas === undefined) return;

    const mol = molCanvas.getMolecule();
    const div = document.getElementById('moleculeImage');
    const width = Math.round(div.offsetWidth);
    const height = Math.round(div.offsetHeight);
    const maxBondLen = getMaxBondLen(mol, width, height);
    mol.scaleToAverageBondLength(maxBondLen);
    molCanvas.resize(width,height);
    molCanvas.loadMolecule(mol);
    paintMolecule(mol);
}

function getMaxBondLen(mol, width, height){
    const dim = mol.getDimension();
    const padding = 10;

    const maxWidth = width - padding;
    const maxHeight = height - padding;

    const scale = Math.min(
        maxWidth / dim.x,
        maxHeight / dim.y
    );

    const safescale = scale * 0.9;

    const current = mol.getAverageBondLength();

    const target = current * safescale;

    return target;
}

function paintMolecule(mol){
        const colorCanvas = document.getElementById('moleculeImageCanvas');
        const ctx = colorCanvas.getContext("2d");
        ctx.globalCompositeOperation='destination-over'; 
        for(const bond of mol.bonds){
            if(!bond.tmpid) bond.tmpid = 0;                                         // could be buggy !!
            if(!bond.a1.tmpid) bond.a1.tmpid = 0;   
            if(!bond.a2.tmpid) bond.a2.tmpid = 0;   
            if(!mol.highlightBonds.includes(bond.tmpid)) continue;
            bond.a1.backgroundColor = mol.highlightAtomColors[bond.a1.tmpid];
            bond.a2.backgroundColor = mol.highlightAtomColors[bond.a2.tmpid];
            bond.backgroundColor = mol.highlightAtomColors[bond.a1.tmpid];
            colorAtom(bond.a1);
            colorAtom(bond.a2);
            colorBond(bond);
        }
        ctx.fillStyle = "white"; 
        ctx.fillRect(0, 0, colorCanvas.width, colorCanvas.height);
}   

function colorAtom(atom, saturation = 0){ 
    const ctx = document.getElementById('moleculeImageCanvas').getContext("2d");  
    let label = atom.label;
    if (!(atom.altLabel ===  undefined)) {
        label = atom.altLabel;
    }
    if (label === 'C') return; 
    const color = atom.backgroundColor;
    const x = atom.x;
    const y = atom.y;
    
    //if (atom.getImplicitHydrogenCount() > 0) label = `${atom.label}H2`;
    const h = molCanvas.styles.atoms_font_size_2D;
    ctx.font = `${h}px Arial`;
    const l = ctx.measureText(label).width;
    
    ctx.strokeStyle = `rgb(${color[0] * 255 + saturation}, ${color[1] * 255 + saturation}, ${color[2] * 255 + saturation})`; 
    ctx.lineWidth = h;  
    ctx.lineCap = "round";
    ctx.beginPath();
    ctx.moveTo(x - l/2 + h * 0.25, y);
    ctx.lineTo(x + l/2 - h * 0.25, y);
    ctx.stroke();
}

function colorBond(bond, saturation = 0){  
    const ctx = document.getElementById('moleculeImageCanvas').getContext("2d");
    const sX = bond.a1.x;
    const sY = bond.a1.y;
    const eX = bond.a2.x;
    const eY = bond.a2.y;
    const color = bond.backgroundColor;
    ctx.strokeStyle = `rgb(${color[0] * 255 + saturation}, ${color[1] * 255 + saturation}, ${color[2] * 255 + saturation})`; 
    ctx.lineWidth = molCanvas.styles.bonds_width_2D * 8;  
    ctx.lineCap = "round";
    ctx.beginPath();
    ctx.moveTo(sX, sY);
    ctx.lineTo(eX, eY);
    ctx.stroke();
}  

function hideLabels(){
    const mol = molCanvas.getMolecule();
    // if atom has a label, set it to undefined such that it is not displayed anymore
    // if atom has a undefined label, set it back to its orig label such that it is displayed again 
    mol.atoms.map(a => a.altLabel ? a.altLabel = undefined : 
                    ( a.altLabel === undefined ? a.altLabel = a.residue : a.altLabel));

    molCanvas.repaint(mol);
    if(mol.selected){
        selectMol(mol.selected, mol.name);
    } else {
        paintMolecule(mol); 
    }
}

function selectMol(nid){
    const mol = molCanvas.getMolecule();
    mol.selected = nid;
    
    // find atoms of selected monomer node in molecule 
    const nodeLabels = [];
    nid.forEach(n => nodeLabels.push(nodeIdLabel.get(n)));
    const selected_atoms = new Set();
    for(const label of nodeLabels){
        mol.monomers[label].forEach(a => selected_atoms.add(a));
    }
    
    // paint molecule in pale color and selected residue in standard color
    let paleStyle = new ChemDoodle.structures.Styles();
    paleStyle.atoms_color = '#D3D3D3';
    paleStyle.bonds_color = '#D3D3D3';
    paleStyle.bonds_width_2D = .9;
    paleStyle.bonds_saturationWidthAbs_2D = 2.6;
    paleStyle.bonds_hashSpacing_2D = 2.5;
    paleStyle.atoms_font_size_2D = 10;
    let standardStyle = new ChemDoodle.structures.Styles();
    standardStyle.atoms_color = 'black';
    standardStyle.bonds_color = 'black';
    standardStyle.bonds_width_2D = .9;
    standardStyle.bonds_saturationWidthAbs_2D = 2.6;
    standardStyle.bonds_hashSpacing_2D = 2.5;
    standardStyle.atoms_font_size_2D = 10;

    let visited = new Set();
    const selected_bonds = new Set();

    for(const bond of mol.bonds){
        const nonColorBond = !mol.highlightBonds.includes(bond.tmpid)
        if(nonColorBond){
            bond.styles = paleStyle;
            continue;
        }

        if(selected_atoms.has(Number(bond.a1.tmpid)) && selected_atoms.has(Number(bond.a2.tmpid))){
            bond.a1.styles = standardStyle;
            bond.a2.styles = standardStyle;
            bond.styles = standardStyle;
            visited.add(bond.a1);
            visited.add(bond.a2);
            selected_bonds.add(bond);
        } else {
            if(!visited.has(bond.a1)) bond.a1.styles = paleStyle;
            if(!visited.has(bond.a2)) bond.a2.styles = paleStyle;
            bond.styles = paleStyle;
        }
    }

    molCanvas.repaint(mol);
    // selected bonds must be painted first to be in the first layer
    for(const bond of selected_bonds){
        colorAtom(bond.a1);
        colorAtom(bond.a2);
        colorBond(bond);
    }

    for(const bond of mol.bonds){
        if(!selected_atoms.has(Number(bond.a1.tmpid)) && !selected_atoms.has(Number(bond.a2.tmpid))){
            if(!mol.highlightBonds.includes(bond.tmpid)) continue;
            colorAtom(bond.a1, 60);
            colorAtom(bond.a2, 60);
            colorBond(bond, 60);
        }
    }
    // set background color back to white
    const colorCanvas = document.getElementById('moleculeImageCanvas');
    const ctx = colorCanvas.getContext("2d");
    ctx.fillStyle = "white"; 
    ctx.fillRect(0, 0, colorCanvas.width, colorCanvas.height);
}

function deselectMol(){
    const mol = molCanvas.getMolecule();
    mol.selected = false;
    let st = new ChemDoodle.structures.Styles();
    st.atoms_color = 'black';
    st.bonds_color = 'black';
    st.bonds_width_2D = .9;
    st.bonds_saturationWidthAbs_2D = 2.6;
    st.bonds_hashSpacing_2D = 2.5;
    st.atoms_font_size_2D = 10;
    
    for(const bond of mol.bonds){
        bond.a1.styles = st;
        bond.a2.styles = st;
        bond.styles = st;
    }
    molCanvas.repaint(mol);
    paintMolecule(mol);
}
function enableZoomMol(el,mol){
    // zoom functionality
    el.addEventListener('wheel', e => e.preventDefault(), {passive: false}); 
    let scalefactor = 1;
    let zoom = e => {
        scalefactor += e.deltaY * -0.01;
        scalefactor = Math.max(scalefactor, 1);
        scalefactor = Math.min(2, scalefactor);
        el.style.transition = "transform 0.1s";
        el.style.transform = `scale(${scalefactor})`;
    }
    el.onwheel = zoom;

    // move functionality
    let dragging = false;
    let offsetX = 0;
    let offsetY = 0;
    let startX, startY;

    el.addEventListener('mousedown', e => {
        dragging = true;
        startX = e.clientX;
        startY = e.clientY;
    }
    );
    el.addEventListener('mousemove', e => {
        if(dragging && scalefactor != 1){
            const dx = e.clientX - startX;
            const dy = e.clientY - startY;

            offsetX += dx;
            offsetY += dy;
            el.style.transform = `translate(${offsetX}px, ${offsetY}px) scale(${scalefactor})`;

            startX = e.clientX;
            startY = e.clientY;
        }
    });
    el.addEventListener('mouseup', e => 
        dragging = false
    );
    el.addEventListener('mouseleave', e => 
        dragging = false
    );

}
// -- deselect handling --
function deselect(){
    if(network === undefined) return;

    // deselect module viewer
    const svgAlign = document.getElementById("moduleViewerAlign");
    const svgBGC = document.getElementById("moduleViewerBGC");
    for(const el of svgBGC.childNodes){
        if(!el.id || !el.displayed) continue;
        el.setAttribute("opacity", "1");
    }
    for(const el of svgAlign.childNodes){
        if(!el.id || !el.displayed) continue;
        el.setAttribute("opacity", "1");
    }

     // deselect alignment Table
    const alignmentTableBody = document.querySelector('#alignmentTable tbody');
    [...alignmentTableBody['rows']].forEach(r => r.style.cssText = "")

    // deselect graph
    network.selectNodes([]); 
    increaseTranparency([]);

    // deselect molecule
    deselectMol();
    
    // deselect nerpa MS if in nerpa MS report 
    if(document.title === 'Nerpa MS Report'){
        deselectMS();
    }
}
// -- select handling --
function select(nid){
    if(network === undefined) return;

    const svgAlign = document.getElementById("moduleViewerAlign");
    const svgBGC = document.getElementById("moduleViewerBGC");
    let nidString = nid.map(String);
    let nidInt = nid.map(Number);
    const hasInvalidIdx = nidInt.some(Number.isNaN);

    // select module viewer
    const cM = Array.from(svgBGC.childNodes).find(m => m.getAttribute('clicked') === 'true');
    const moduleClicked = cM != undefined;
    const clickedGene = moduleClicked ? splitId(cM.id).gid : '';
    const noNid = nid.includes('---');

    for(const el of svgBGC.childNodes){
        if(!el.id) continue;
        const split = splitId(el.id);
        const noNidCondition = moduleClicked && noNid && split.gid != clickedGene;
        if(split.nid.some(n => nidString.includes(n)) && !noNidCondition){
            el.setAttribute("opacity", "1");
        } else {
            el.setAttribute("opacity", "0.3");
        }
    }
    for(const el of svgAlign.childNodes){
        if(!el.id) continue;
        const split = splitId(el.id);
        const noNidCondition = moduleClicked && noNid && split.gid != clickedGene;
        if(split.nid.some(n => nidString.includes(n)) && !noNidCondition){
            el.setAttribute("opacity", "1");
        } else {
            el.setAttribute("opacity", "0.3");
        }
    }
    // select alignment Table
    const alignmentTableBody = document.querySelector('#alignmentTable tbody');
    [...alignmentTableBody['rows']].forEach(r => {
        const rid = r.cells[6].textContent;
        if (nidString.includes(rid)){
            r.style.border = "2px solid #9EC37B";
            r.style.background = '#ddfcdb';
        } else {
            r.style.cssText = "";
        }
    });

    if(hasInvalidIdx) {
        nidString = nidInt.filter(nid => !Number.isNaN(nid)).map(String);
    } 
    if(nidString.length === 0){
        // select graph
        network.selectNodes([]);  
        increaseTranparency(nidInt);
    } else {
        // select graph
        network.selectNodes(nidString);  
        increaseTranparency(nidString);
    }
    // select molecule only if nerpa Report, in nerp MS no connection betwenn graph and molecule
    if(document.title === 'Nerpa Report'){
        selectMol(nidString);
    }
    
}
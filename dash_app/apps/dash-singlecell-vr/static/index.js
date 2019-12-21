/* eslint-disable no-unused-vars */
/*global Fuse, THREE, AFRAME, JSZip*/

let report = {};
let freeMove = true;
let positionIndex = 0;
let currentBranch = null;
let cameraTrajectory = null;
const branchClasses = [];
const cellColors = {};

document.getElementById("moveToggle").addEventListener("click", () => {
  freeMove = !freeMove;
});

const unzip = async (uuid) => {
  const zipper = new JSZip();
  const response = await fetch('/download/' + uuid + '.zip');
  const blob = await response.blob();
  const result = await zipper.loadAsync(blob)
  return result;
}

const mobilecheck = () => {
  var check = false;
  // eslint-disable-next-line no-useless-escape
  (function(a){if(/(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino/i.test(a)||/1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas\-|your|zeto|zte\-/i.test(a.substr(0,4))) check = true;})(navigator.userAgent||navigator.vendor||window.opera);
  return check;
};

const visibleHeightAtZDepth = ( depth ) => {
  const camera = AFRAME.scenes[0].camera;
  // compensate for cameras not positioned at z=0
  const cameraOffset = camera.position.z;
  if ( depth < cameraOffset ) depth -= cameraOffset;
  else depth += cameraOffset;

  // vertical fov in radians
  const vFOV = camera.fov * Math.PI / 180;

  // Math.abs to ensure the result is always positive
  return 2 * Math.tan( vFOV / 2 ) * Math.abs( depth );
};

const visibleWidthAtZDepth = ( depth ) => {
  const camera = AFRAME.scenes[0].camera;
  const height = visibleHeightAtZDepth( depth, camera );
  let width = height * camera.aspect;
  return width;
};

const setHudPosition = ( fovWidth, fovHeight, depth) => {
  document.getElementById('hud').object3D.position.set(-fovWidth/2 + .25, fovHeight/2 - .25, depth);
}

setHudPosition(visibleWidthAtZDepth(-1), visibleHeightAtZDepth(-1), -1);

const summonMenu = () => {
  const camera = document.getElementById("curve-camera");
  const start = new THREE.Vector3();
  camera.object3D.getWorldPosition(start);
  const direction = new THREE.Vector3(0, 0, -1);
  direction.applyQuaternion(camera.object3D.quaternion);
  const newPos = new THREE.Vector3();
  newPos.addVectors( start, direction.multiplyScalar( 5 ) );
  const menu = document.getElementById("menu").object3D;
  menu.position.set(newPos.x, newPos.y, newPos.z);
}

const initializeGui = () => {
  document.getElementById("test_input").setAttribute('value', "");
  document.getElementById("result1").setAttribute('value', "");
  document.getElementById("result2").setAttribute('value', "");
  document.getElementById("result3").setAttribute('value', "");
}
initializeGui();

let geneList = [];

const movement = (num) => {
  let direction = new THREE.Vector3();
  const camera = AFRAME.scenes[0].camera;
  camera.getWorldDirection( direction );
  direction.multiplyScalar(num);
  const cameraEl = document.getElementById('rig');
  var pos = cameraEl.getAttribute("position");
  pos.x += direction.x
  pos.y += direction.y
  pos.z += direction.z
  const mapPlayer = document.getElementById('mapPlayer').object3D;
  mapPlayer.position.set((pos.x + direction.x)  * .01, (pos.y + direction.y) * .01, (pos.z + direction.z) * .01);
}

let currentSearch = '';
document.body.addEventListener('keydown', (e) => {
  var options = {
    shouldSort: true,
    threshold: 0.6,
    location: 0,
    distance: 100,
    maxPatternLength: 32,
    minMatchCharLength: 1,
    keys: [
      "gene",
    ]
  };
  const resultsEntity = document.getElementById("test_input");
  let result1 = '';
  let result2 = '';
  let result3 = '';
  if (e.code === 'Space') {
    summonMenu();
    currentSearch = '';
  } else if (e.key === "Shift") {
    if (mobilecheck()) {
      const hud = document.getElementById("hud").object3D;
      if (!hud.visible) {
        hud.position.set(0, 0, -.5);
        hud.visible = true;
      }
      document.getElementById("cursor").object3D.visible = false;
    }
  } else if (e.keyCode === 38) {
    if (freeMove) {
      movement(.05);
    } else {
      moveCamera(e);
    }
  } else if (e.keyCode === 40) {
    if (freeMove) {
      movement(-.05);
    } else {
      moveCamera(e);
    }
  } else if (e.code === 'Enter') {
    currentSearch = '';
    viewGene('metadata', 'label_color');
  } else if (e.key.length === 1) {
    currentSearch = currentSearch + e.key;
    let fuse = new Fuse(geneList, options);
    let result = fuse.search(currentSearch);
    result1 = result.length > 0 ? result[0].gene : "";
    result2 = result.length > 1 ? result[1].gene : "";
    result3 = result.length > 2 ? result[2].gene : "";

  }
  resultsEntity.setAttribute('text', 'value', "Search: " + currentSearch);
  document.getElementById("result1").setAttribute('text', 'value', result1);
  document.getElementById("result2").setAttribute('text', 'value', result2);
  document.getElementById("result3").setAttribute('text', 'value', result3);
});

document.body.addEventListener('keyup', (e) => {
  if (e.key === "Shift") {
    if (mobilecheck()) {
      const hud = document.getElementById("hud").object3D;
      if (hud.visible) {
        hud.visible = false;
      }
      document.getElementById("cursor").object3D.visible = true;
    }
  }
});

const resultElements = ["result1", "result2", "result3"];
resultElements.forEach((element) => {
  const result = document.getElementById(element);
  result.addEventListener("click", () => {
    viewGene('gene_' + result.getAttribute('text').value, 'color');
  });
});

setInterval(() => {
    const drawContainer = document.getElementById("drawContainer");
    const drawContainerRotation = drawContainer.object3D.rotation;
    const hudMapContainer = document.getElementById("hudMapContainer");
    hudMapContainer.object3D.rotation.set(drawContainerRotation._x, drawContainerRotation._y, drawContainerRotation._z);
}, 50);

document.getElementById("pauseGlobalRotation").addEventListener("click", () => {
  const drawContainer = document.getElementById("drawContainer");
  const isRotating = drawContainer.isPlaying;
  if (isRotating) {
    drawContainer.pause();
    drawContainer.setAttribute("rotation", "0 0 0");
  } else {
    drawContainer.play();
  }
});

const viewGene = async (geneFileName, colorField) => {
  const gene = await report.file(geneFileName + ".json").async("string");
  const cellsByGene = JSON.parse(gene);
  const cellsContainer = document.getElementById("cells");
  cellsByGene.forEach((cell) => {
    const cellElement = document.getElementById(cell.cell_id);
    cellElement.setAttribute("color", cell[colorField]);
  });
}

document.querySelector('a-scene').addEventListener('enter-vr', () => {
  const hud = document.getElementById("hud");
  hud.setAttribute('material', 'color', 'white');
  const labels = document.getElementById("curve-labels").childNodes;
  labels.forEach((label) => {
    label.setAttribute("text", "color", "white");
  })
  setHudPosition(visibleWidthAtZDepth(-1) - .5, visibleHeightAtZDepth(-1), -1);
  if (mobilecheck()) {
    document.getElementById('hud').object3D.visible = false;
  }
  const legend = document.getElementById('legend');
  if (legend !== null) {
    document.getElementById('legend').setAttribute('panel-color', 'black');
  }
});


document.querySelector('a-scene').addEventListener('exit-vr', () => {
  const hud = document.getElementById("hud");
  hud.setAttribute('material', 'color', 'gray');
  const labels = document.getElementById("curve-labels").childNodes;
  labels.forEach((label) => {
    label.setAttribute("text", "color", "black");
  })
  setHudPosition(visibleWidthAtZDepth(-1), visibleHeightAtZDepth(-1), -1);
  const legend = document.getElementById('legend');
  if (legend !== null) {
    document.getElementById('legend').setAttribute('panel-color', 'white');
  }
});

const getZMax = (curves) => {
    let maxZ = Number.NEGATIVE_INFINITY;
    curves.forEach((curve) => {
      curve.xyz.forEach((coord) => {
        if (typeof coord.z !== 'undefined' && Math.abs(coord.z) > maxZ) {
            maxZ = Math.abs(coord.z);
        }
      });
    });
    return maxZ * 100;
}

const getMedian = (values) => {
    const sorted = [...values].sort();
    return sorted[Math.floor(sorted.length/2)];
}

const getCameraTrajectory = (branches) => {
  currentBranch = branches[0].branch_id;
  const branchPositions = {};
  branches.forEach((branch) => {
    const positions = [];
    branch.xyz.forEach((coord, _) => {
        const position = `${coord.x} ${coord.y} ${coord.z}`;
        positions.push(position);
    });
    branchPositions[branch.branch_id] = positions;
    branchClasses.push('.' + branch.branch_id);
  });
  return branchPositions;
};

const groupBy = (list, keyGetter) => {
    const map = new Map();
    list.forEach((item) => {
        const key = keyGetter(item);
        const collection = map.get(key);
        if (!collection) {
            map.set(key, [item]);
        } else {
            collection.push(item);
        }
    });
    return map;
}

const getFileText = (name) => {
  return fetch(name + '.json')
    .then(response => response.text())
    .then(text => {
        return JSON.parse(text);
    });
}

const createCellMetadataObject = (metadata) => {
  const cellObjects = {};
  metadata.forEach((cell) => {
    cellObjects[cell.cell_id] = {"label": cell.label, "label_color": cell.label_color, "cluster_color": cell.cluster_color}
  });
  return cellObjects;
}

const renderPagaCells = (cells, cellMetadata) => {
  const cellEntities = Array.from(cells.map((cell) => {
    const x = cell.x * .0004;
    const y = cell.y * .0004;
    const color = cellMetadata[cell.cell_id].cluster_color;
    return `<a-sphere id="${cell.cell_id}" position="${x} ${y} -1" radius=".03" color="${color}"></a-sphere>`
  }));
  document.getElementById('pagacells').innerHTML = cellEntities.join(" ");
}

const setInitialCameraPositionPaga = (nodes) => {
  const yValues = Array.from(Object.values(nodes).map(node => node.xy.y * .0004));
  const xValues = Array.from(Object.values(nodes).map(node => node.xy.x * .0004));
  const xMax = Math.max(...xValues);
  const xMin = Math.min(...xValues);
  const xRange = xMax - xMin;
  const yMax = Math.max(...yValues);
  const yMin = Math.min(...yValues);
  const xMidpoint = (xMax + xMin) / 2;
  const yMidpoint = (yMax + yMin) / 2;

  const camera_el = document.getElementById("rig");
  camera_el.object3D.position.set(xMidpoint, yMidpoint, xRange + 1);
}

const renderLegend = (metadata) => {
  const legendColors = {};
  metadata.forEach((metadatum) => {
    legendColors[metadatum.cluster] = metadatum.cluster_color;
  });
  const legend = document.getElementById('legend');
  Object.keys(legendColors).forEach((key) => {
    const el = document.createElement("a-gui-button");
    el.setAttribute("width", "2.5");
    el.setAttribute("height", ".25");
    el.setAttribute("value", key);
    el.setAttribute("font-color", "black");
    el.setAttribute("background-color", legendColors[key]);
    legend.appendChild(el);
  });
}

const renderPaga = (edges, nodes, scatter, metadata) => {
  setInitialCameraPositionPaga(nodes);
  renderLegend(metadata);
  const branches = [];
  const edgeWeights = {};
  edges.forEach((edge, _) => {
      const edgeId = edge.nodes[0] + '_' + edge.nodes[1];
      if (!branches.includes(edgeId)) {
          branches.push(edgeId);
      }
    edgeWeights[edgeId] = edge.weight;
  });
  const [branch_els, branch_draw_els] = createCurveEnities(branches);
  // setDrawContainerContent(branch_els, branch_draw_els);
  const clusterColors = createCellMetadataObject(metadata);
  const cell_el = document.getElementById("cells");
  const cellEntities = [];
  const nodePositions = {};
  Object.values(nodes).forEach((cell_point, _) => {
    let x = cell_point.xy.x * .0004;
    let y = cell_point.xy.y * .0004;
    const stream_cell = `<a-sphere text="value: ${cell_point.node_name}; width: 6; color: black; align: center; side: double; zOffset: .1" id="${cell_point.node_name}" position="${x} ${y} -1" color="${clusterColors[cell_point.node_name]}" radius=".1" billboard></a-sphere>`;
    cellEntities.push(stream_cell);
    nodePositions[cell_point.node_name.replace(/\D/g,'')] = {"x": x, "y": y, "z": -1};
  });
  cell_el.innerHTML = cellEntities.join(" ");
  const thickLines = [];
  branch_els.forEach((branch) => {
    const curveref = branch.split(" ");
    // TODO: Bad, evil, fix!
    const curveid = curveref[1].split("=")[1];
    const [startNode, endNode] = strip(curveid).split("_");
    const thickLine = `<a-entity meshline="lineWidth: ${edgeWeights[strip(curveid)] * 10}; path: ${nodePositions[startNode].x} ${nodePositions[startNode].y} ${nodePositions[startNode].z}, ${nodePositions[endNode].x} ${nodePositions[endNode].y} ${nodePositions[endNode].z}; color: black"></a-entity>`
    thickLines.push(thickLine);
  });
  document.getElementById("thicklines").innerHTML = thickLines.join(" ");
  document.getElementById("thicklinesMap").innerHTML = thickLines.join(" ");
  renderPagaCells(scatter, clusterColors);
}

const strip = (str) => {
    return str.replace(/^\"+|\"+$/g, '');
}

const createBranchPoints = (curve) => {
  const curvePoints = [];
  const midpoint = curve.xyz[Math.floor(curve.xyz.length / 2)];
  const labelEntity = document.createElement("a-entity");
  const textValue = `value: ${curve.branch_id}; color:black; align: center; side: double; width: 6`;
  labelEntity.setAttribute("text", textValue);
  const labelPosition = `${midpoint.x * 100} ${midpoint.y * 100} ${midpoint.z * 100}`;
  labelEntity.setAttribute("position", labelPosition);
  labelEntity.setAttribute("billboard", "");
  labelEntity.className = curve.branch_id;
  labelEntity.addEventListener('click', () => {
    currentBranch = curve.branch_id;
  })
  const curveLabels = document.getElementById('curve-labels');
  curveLabels.appendChild(labelEntity);
  curve.xyz.forEach((coord, _) => {
      const curvePoint = `<a-curve-point position="${coord.x * 100} ${coord.y * 100} ${coord.z * 100}"></a-curve-point>`;
      curvePoints.push(curvePoint);
  });
  return curvePoints;
}

const setInitialCameraPosition = (curves) => {
  const zMax = getZMax(curves);
  const yValues = Array.from(curves.flatMap((curve) => {
    return Array.from(curve.xyz.map(coord => coord.y));
  }));
  const xValues = Array.from(curves.flatMap((curve) => {
    return Array.from(curve.xyz.map(coord => coord.x))
  }));
  const yMedian = getMedian(yValues) * 100;
  const xMedian = getMedian(xValues) * 100;

  const camera_el = document.getElementById("rig");
  camera_el.object3D.position.set(xMedian, yMedian, zMax + 1.2);
  const mapPlayer = document.getElementById("mapPlayer");
  mapPlayer.object3D.position.set(xMedian * .01, yMedian * .01, (zMax + 1.2) * .01)
}

const createCurveEnities = (branches) => {
  const branch_els = [];
  const branch_draw_els = [];
  branches.forEach((branch, _) => {
      const branch_el = `<a-curve id="${branch}" ></a-curve>`;
      branch_els.push(branch_el);
      const branch_draw_el = `<a-draw-curve cursor-listener curveref="#${branch}" material="shader: line; color: blue;" geometry="primitive: " ></a-draw-curve>`;
      branch_draw_els.push(branch_draw_el);
  });
  return [branch_els, branch_draw_els];
}

const setDrawContainerContent = (branch_els, branch_draw_els) => {
  const branch_container_el = document.getElementById("curve-container");
  branch_container_el.innerHTML = branch_els.join(" ");
  const map_branch_container = document.getElementById("curve-map");
  map_branch_container.innerHTML = branch_els.join(" ");
  const branch_draw_container = document.getElementById("curve-draw");
  branch_draw_container.innerHTML = branch_draw_els.join(" ");
  const map_draw_container = document.getElementById("draw-map");
  map_draw_container.innerHTML = branch_draw_els.join(" ");
}

const renderStream = async (curves, cells, metadata) => {
  setInitialCameraPosition(curves);
  const camvec = new THREE.Vector3();
  const camera = AFRAME.scenes[0].camera;
  camera.getWorldPosition(camvec);
  document.getElementById('draw-map').object3D.lookAt(-camera.position.x, -camera.position.y, -camera.position.z);
  const branches = [];
  curves.forEach((coord, _) => {
      if (!branches.includes(coord.branch_id)) {
          branches.push(coord.branch_id);
      }
  });

  const [branch_els, branch_draw_els] = createCurveEnities(branches);

  setDrawContainerContent(branch_els, branch_draw_els);

  curves.forEach((curve) => {
    const points = createBranchPoints(curve);
    const branch_el = document.getElementById(curve.branch_id);
    branch_el.innerHTML = points.join(" ")
  })

  metadata.forEach((cell) => {
    cellColors[cell.cell_id] = cell.label_color;
  });

  renderStreamCells(cells);
}

const renderStreamCells = async (cells) => {
  const cell_el = document.getElementById("cells");
  const cellEntities = [];
  cells.forEach((cell_point, _) => {
    const stream_cell = `<a-sphere id="${cell_point.cell_id}" position="${cell_point.x * 100} ${cell_point.y * 100} ${cell_point.z * 100}" color="${cellColors[cell_point.cell_id]}" radius=".05" shadow></a-sphere>`;
    cellEntities.push(stream_cell);
  });
  cell_el.innerHTML = cellEntities.join(" ");
}

let clickCount = 1;

const move = (position) => {
    const camera_el = document.getElementById("rig");
    let positionSplit = position.split(" ");
    const scaled = positionSplit.map((coord) => {
        return coord * 100;
    });
    camera_el.object3D.position.set(...scaled);
}


const moveCamera = (e) => {
  if (e.keyCode === 38) {
    if (positionIndex !== cameraTrajectory[currentBranch].length) {
      move(cameraTrajectory[currentBranch][positionIndex]);
      positionIndex += 1;
    }
  } else if (e.keyCode === 40) {
    if (positionIndex >= 0) {
      positionIndex = Math.max(positionIndex-1, 0);
      move(cameraTrajectory[currentBranch][positionIndex])
    }
  }
}

// TODO: Quick and dirty. strips everything besides gui-interactable. Fix Later
const makeIntersectable = (objects) => {
  const cursor = document.getElementById("cursor");
  const currentIntersectable = cursor.getAttribute('raycaster');
  cursor.setAttribute('raycaster', 'objects', '[gui-interactable], ' + objects.join(", "));
  console.log(cursor.getAttribute('raycaster').objects);
}

const getGeneList = (report) => {
  const allFileNames = Object.keys(report.files);
  const geneNames = [];
  allFileNames.forEach((file) => {
    const splitName = file.split("_");
    if (splitName.length > 1 && splitName[0] === 'gene') {
      geneNames.push({'gene': splitName[1].split('.')[0]})
    }
  });
  return geneNames;
}

const initialize = async (uuid) => {
  const result = await unzip(uuid);
  report = result;
  if (Object.keys(result.files).includes("paga_nodes.json")) {
    const edges = await result.file("paga_edges.json").async("string");
    const nodes = await result.file("paga_nodes.json").async("string");
    const scatter = await result.file("scatter.json").async("string");
    const metadata = await result.file("metadata.json").async("string");
    document.getElementById("moveToggle").remove();
    renderPaga(JSON.parse(edges), JSON.parse(nodes), JSON.parse(scatter), JSON.parse(metadata));
  } else {
    const streamFile = await result.file("stream.json").async("string");
    const scatterFile = await result.file("scatter.json").async("string");
    const metadataFile = await result.file("metadata.json").async("string");
    cameraTrajectory = getCameraTrajectory(JSON.parse(streamFile));
    makeIntersectable(branchClasses);
    document.getElementById('legend').remove();
    document.getElementById("moveToggle").remove();
    renderStream(JSON.parse(streamFile), JSON.parse(scatterFile), JSON.parse(metadataFile));
  }
  geneList = getGeneList(result);
}

window.onload = () => {
  const uuid = window.location.href.split("/").slice(-1)[0];
  initialize(uuid);
}

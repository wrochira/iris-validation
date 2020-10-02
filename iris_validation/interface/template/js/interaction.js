/*
Copyright 2020 William Rochira at York Structural Biology Laboratory
*/

const COLORS = {  'L_GREY' : 'rgb(150, 150, 150)',
                  'VL_GREY' : 'rgb(200, 200, 200)',
                  'WHITE' : 'rgb(255, 255, 255)',
                  'BAR_GREEN' : 'rgb(90, 237, 141)',
                  'BAR_ORANGE' : 'rgb(247, 212, 134)',
                  'BAR_RED' : 'rgb(240, 106, 111)' };


function mean(values) {
  let sum = values.reduce(function(sum, value) {
    return sum + value;
  }, 0);
  let avg = sum / values.length;
  return avg;
};


function standardDeviation(values) {
  avg = mean(values);
  let squareDiffs = values.map(function(value) {
    var diff = value - avg;
    var sqrDiff = diff * diff;
    return sqrDiff;
  });
  let variance = mean(squareDiffs);
  let stdDev = Math.sqrt(variance);
  return stdDev;
};


function percentile(values, p) {
  values = values.sort();
  var pos = ((values.length) - 1) * p;
  var base = Math.floor(pos);
  var rest = pos - base;
  if( (values[base+1]!==undefined) ) {
    return values[base] + rest * (values[base+1] - values[base]);
  } else {
    return values[base];
  }
};


function coordsFromAngle(centre, angle, pointRadius) {
  let xc = centre[0];
  let yc = centre[1];
  let x1 = xc + pointRadius * Math.sin(angle);
  let y1 = yc - pointRadius * Math.cos(angle);
  return [Math.round(x1, 2), Math.round(y1, 2)];
};


function scrollToId(elementId) {
  const el = document.getElementById(elementId);
  const offset = 60;
  const bodyRect = document.body.getBoundingClientRect().top;
  const elRect = el.getBoundingClientRect().top;
  const elPos = elRect - bodyRect;
  const offsetPos = elPos - offset;

  window.scrollTo({
    top: offsetPos,
    behavior: 'smooth'
  });
};


function toggleModel() {
  if (document.getElementById('model-toggle').checked) {
    selectedModel = 1;
  } else {
    selectedModel = 0;
  }
  // If selected residue is null on the newly-selected model, cycle residues
  while (discreteMetrics[selectedModel][selectedChain][selectedResidue] === null) {
    selectedResidue++;
    selectedResidue = selectedResidue % chainLengths[selectedChain];
    setSelector(selectedChain, selectedResidue);
  };
  setBinarySegments(selectedModel);
  setClashMarkers(selectedModel);
  setShade(selectedModel);
  animatePolyline(selectedModel);
  //setRadarChart(selectedModel, selectedChain, selectedResidue);
  setResidueChart(selectedModel, selectedChain, selectedResidue);
  setResidueChartRanges();
};


function setChain(chainID) {
  isDragging = false;
  selectedChain = chainID;
  selectedResidue = 0;
  //clearRadarChart();
  //clearResidueChart();
  setResidueChart(selectedModel, selectedChain, selectedResidue);
  setSelector(selectedChain, selectedResidue);
  setIrisChart(selectedChain);
};


function setIrisChart(chainID) {
  for (var i = 0; i < numChains; ++i) {
    if (i === chainID) {
      document.getElementById('iris-chart-' + i).style.display = '';
      document.getElementById('chain-button-' + i).style['color'] = COLORS['L_GREY'];
    } else {
      document.getElementById('iris-chart-' + i).style.display = 'none';
      document.getElementById('chain-button-' + i).style['color'] = null;
    };
  };
};


/*
function setChartComponents(mode) {
  let components = [ 'discrete', 'markers', 'shade' ];
  for (var componentID = 0; componentID < components.length; ++componentID) {
    for (var chainID = 0; chainID < numChains; ++chainID) {
      let chartID = 'iris-chart-' + chainID;
      let latest = document.getElementById(chartID + '-' + components[componentID] + '-1');
      let previous = document.getElementById(chartID + '-' + components[componentID] + '-0');
      if (latest != null && previous != null) {
        if (mode === 0) {
          latest.setAttribute('opacity', 0);
          previous.setAttribute('opacity', 1);
        } else if (mode === 1) {
          latest.setAttribute('opacity', 1);
          previous.setAttribute('opacity', 0);
        };
      };
    };
  };
};
*/


function setBinarySegments(mode) {
  for (var chainID = 0; chainID < numChains; ++chainID) {
    let chartID = 'iris-chart-' + chainID;
    for (var metricID = 0; metricID < numMetrics; ++metricID) {
      let latestBSG = document.getElementById(chartID + '-discrete-' + metricID + '-1');
      let previousBSG = document.getElementById(chartID + '-discrete-' + metricID + '-0');
      if (latestBSG != null && previousBSG != null) {
        if (mode === 0) {
          latestBSG.setAttribute('opacity', 0);
          previousBSG.setAttribute('opacity', 1);
        } else if (mode === 1) {
          latestBSG.setAttribute('opacity', 1);
          previousBSG.setAttribute('opacity', 0);
        };
      }
    };
  };
};


function setClashMarkers(mode) {
  for (var chainID = 0; chainID < numChains; ++chainID) {
    let chartID = 'iris-chart-' + chainID;
    let latestMarkers = document.getElementById(chartID + '-markers-1');
    let previousMarkers = document.getElementById(chartID + '-markers-0');
    if (latestMarkers != null && previousMarkers != null) {
      if (mode === 0) {
        latestMarkers.setAttribute('opacity', 0);
        previousMarkers.setAttribute('opacity', 1);
      } else if (mode === 1) {
        latestMarkers.setAttribute('opacity', 1);
        previousMarkers.setAttribute('opacity', 0);
      };
    };
  };
};


function setShade(mode) {
  for (var chainID = 0; chainID < numChains; ++chainID) {
    let chartID = 'iris-chart-' + chainID;
    let latestShade = document.getElementById(chartID + '-shade-1');
    let previousShade = document.getElementById(chartID + '-shade-0');
    if (latestShade != null && previousShade != null) {
      if (mode === 0) {
        latestShade.setAttribute('opacity', 0);
        previousShade.setAttribute('opacity', 1);
      } else if (mode === 1) {
        latestShade.setAttribute('opacity', 1);
        previousShade.setAttribute('opacity', 0);
      };
    };
  };
};


function animatePolyline(mode) {
  for (var chainID = 0; chainID < numChains; ++chainID) {
    let chartID = 'iris-chart-' + chainID;
    for (var metricID = 0; metricID < numMetrics; ++metricID) {
      let toLatest = document.getElementById(chartID + '-animation-' + metricID + '-1');
      let toPrevious = document.getElementById(chartID + '-animation-' + metricID + '-0');
      if (toLatest != null && toPrevious != null) {
        if (mode === 0) {
          toPrevious.beginElement();
        } else if (mode === 1) {
          toLatest.beginElement();
        };
      };
    };
  };
};


function handleSegment(actionID, residueID) {
  if (actionID === 1 || actionID === 2 && isDragging) {
    isDragging = true;
    selectedResidue = residueID;
    //setRadarChart(selectedModel, selectedChain, selectedResidue);
    setResidueChart(selectedModel, selectedChain, selectedResidue);
    setSelector(selectedChain, selectedResidue);
  } else if (actionID === 3) {
    isDragging = false;
  };
};


function setSelector(chainID, residueID) {
  let chartID = 'iris-chart-' + chainID;
  let selector = document.getElementById(chartID + '-selector');
  selector.setAttribute('transform', 'rotate(' + (360-gapDegrees)/chainLengths[chainID]*residueID + ', 500, 500' + ')');
  let interactionSegments = document.querySelectorAll('[id^=' + chartID + '-intseg-]');
  for (var i = 0; i < interactionSegments.length; ++i) {
    if (i === residueID) {
      interactionSegments[i].setAttribute('stroke-opacity', 0.25);
      interactionSegments[i].setAttribute('fill-opacity', 0.25);
    } else {
      interactionSegments[i].setAttribute('stroke-opacity', 0);
      interactionSegments[i].setAttribute('fill-opacity', 0);
    };
  };
};


function clearRadarChart() {
  let radarChart = document.getElementById('radar-chart');
  let canvas = [parseInt(radarChart.getAttribute('viewBox').split(' ')[2]), parseInt(radarChart.getAttribute('viewBox').split(' ')[3])];
  let centre = [canvas[0]/2, canvas[1]/2];
  let polygonPoints = [ ];
  let circles = document.querySelectorAll('[id^=metric-point-]');
  for (var i = 0; i < circles.length; ++i) {
    polygonPoints.push(centre);
    circles[i].setAttribute('cx', centre[0]);
    circles[i].setAttribute('cy', centre[1]);
    circles[i].setAttribute('fill', COLORS['WHITE']);
    circles[i].setAttribute('fill-opacity', 0);
    circles[i].setAttribute('stroke-opacity', 0);
  };
  let polygon = document.getElementById('connecting-polygon');
  polygon.setAttribute('points', polygonPoints);
  polygon.setAttribute('fill-opacity', 0);
  polygon.setAttribute('stroke-opacity', 0);

  let summaryEl = document.getElementById('residue-summary');
  summaryEl.textContent = 'Select a residue...';
};


function setRadarChart(modelID, chainID, residueID) {
  let aaCode, seqNum;
  let radarChart = document.getElementById('radar-chart');
  let canvas = [parseInt(radarChart.getAttribute('viewBox').split(' ')[2]), parseInt(radarChart.getAttribute('viewBox').split(' ')[3])];
  let centre = [canvas[0]/2, canvas[1]/2];
  let axisRadius = Math.min(canvas[0], canvas[1])/2 - Math.min(canvas[0], canvas[1])/20;
  let polygonPoints = [ ];
  let polygon = document.getElementById('connecting-polygon');
  let circles = document.querySelectorAll('[id^=metric-point-]');
  if (absoluteMetrics[modelID][chainID][residueID] === null) {
    aaCode = 'N/A';
    seqNum = 'N/A';
    for (var i = 0; i < circles.length; ++i) {
      polygonPoints.push(centre);
      circles[i].setAttribute('cx', centre[0]);
      circles[i].setAttribute('cy', centre[1]);
      circles[i].setAttribute('fill', COLORS['WHITE']);
      circles[i].setAttribute('fill-opacity', 0);
      circles[i].setAttribute('stroke-opacity', 0);
    };
    polygon.setAttribute('points', polygonPoints);
    polygon.setAttribute('fill-opacity', 0);
    polygon.setAttribute('stroke-opacity', 0);
  } else {
    aaCode = residueCodes[modelID][chainID][residueID];
    seqNum = sequenceNumbers[modelID][chainID][residueID];
    let percentileValues = percentileMetrics[modelID][chainID][residueID];
    let angleDelta = 2*Math.PI / percentileValues.length;
    for (var i = 0; i < circles.length; ++i) {
      let newCentre = coordsFromAngle(centre, i*angleDelta, axisRadius*percentileValues[i]/100);
      polygonPoints.push(newCentre);
      circles[i].setAttribute('cx', newCentre[0]);
      circles[i].setAttribute('cy', newCentre[1]);
      let color = percentileValues[i] === null ? COLORS['VL_GREY'] : percentileValues[i] < 25 ? COLORS['BAR_RED'] : percentileValues[i] < 75 ? COLORS['BAR_ORANGE'] : COLORS['BAR_GREEN'];
      circles[i].setAttribute('fill', color);
      circles[i].setAttribute('fill-opacity', 0.5);
      circles[i].setAttribute('stroke-opacity', 1);
    };
    polygon.setAttribute('points', polygonPoints);
    polygon.setAttribute('fill-opacity', 0.5);
    polygon.setAttribute('stroke-opacity', 1);
  };

  let summaryEl = document.getElementById('residue-summary');
  summaryEl.textContent = 'Chain ' + String.fromCharCode(65 + chainID) + ', ' + 'Residue ' + seqNum + ' (' + aaCode + ')';
};


function handleRadarChart(actionID, metricID) {
  let floatingLabel = document.getElementById('floating-label-group');
  let floatingLabelTextTop = document.getElementById('floating-label-text-top');
  let floatingLabelTextBottom = document.getElementById('floating-label-text-bottom');
  // On mouse over
  if (actionID === 2) {
    let absoluteValue = absoluteMetrics[selectedModel][selectedChain][selectedResidue][metricID];
    let percentileValue = percentileMetrics[selectedModel][selectedChain][selectedResidue][metricID];
    let metricPoint = document.getElementById('metric-point-' + metricID);
    if (absoluteValue === null) {
      floatingLabelTextTop.textContent = 'N/A';
      floatingLabelTextBottom.textContent = 'N/A';
    } else {
      let roundedAbsolute;
      if (absoluteValue < 1) {
        roundedAbsolute = Math.round(absoluteValue*1000)/1000;
      } else if (absoluteValue < 10) {
        roundedAbsolute = Math.round(absoluteValue*100)/100;
      } else if (absoluteValue < 100) {
        roundedAbsolute = Math.round(absoluteValue*10)/10;
      } else {
        roundedAbsolute = Math.round(absoluteValue);
      }
      floatingLabelTextTop.textContent = roundedAbsolute;
      floatingLabelTextBottom.textContent = '(' + percentileValue + '%)';
    }
    floatingLabel.setAttribute('transform', 'translate(' + (parseInt(metricPoint.getAttribute('cx'))) + ', ' + (parseInt(metricPoint.getAttribute('cy'))-30) + ')');
    floatingLabel.setAttribute('opacity', 1);
  // On mouse out
  } else if (actionID === 3) {
    floatingLabelTextTop.textContent = '';
    floatingLabelTextBottom.textContent = '';
    floatingLabel.setAttribute('opacity', 0);
  };
};


function getResidueChartDims() {
  let barChartsContainer = document.getElementById('bar-charts-container');
  let bccPoints = [ ];
  for (var i=0; i<barChartsContainer.points.numberOfItems; ++i) {
    let x = barChartsContainer.points.getItem(i).x;
    let y = barChartsContainer.points.getItem(i).y;
    bccPoints.push([x, y]);
  };
  barOffsetY = bccPoints[2][1];
  barMultiplierY = -(bccPoints[2][1]-bccPoints[0][1]) / 100;
};


function setResidueChart(modelID, chainID, residueID) {
  // Get data for checkboxes
  let box1Color = '';
  let box2Color = '';
  let box1Text = '';
  let box2Text = '';
  let rama = discreteMetrics[modelID][chainID][residueID][0];
  if (rama === null) {
    box1Color = COLORS['VL_GREY'];
    box1Text = 'N/A';
  } else if (rama == '0') {
    box1Color = COLORS['BAR_RED'];
    box1Text = 'Unfavoured';
  } else if (rama == '1') {
    box1Color = COLORS['BAR_ORANGE'];
    box1Text = 'Allowed';
  } else if (rama == '2') {
    box1Color = COLORS['BAR_GREEN'];
    box1Text = 'Favoured';
  };
  let rota = discreteMetrics[modelID][chainID][residueID][1];
  if (rota === null) {
    box2Color = COLORS['VL_GREY'];
    box2Text = 'N/A';
  } else if (rota == '0') {
    box2Color = COLORS['BAR_RED'];
    box2Text = 'Unfavoured';
  } else if (rota == '1') {
    box2Color = COLORS['BAR_ORANGE'];
    box2Text = 'Allowed';
  } else if (rota == '2') {
    box2Color = COLORS['BAR_GREEN'];
    box2Text = 'Favoured';
  };
  // Set checkbox color and text
  document.getElementById('checkbox-1').setAttribute('fill', box1Color);
  document.getElementById('checkbox-2').setAttribute('fill', box2Color);
  document.getElementById('checkbox-1-text').textContent = box1Text;
  document.getElementById('checkbox-2-text').textContent = box2Text;
  // Make box plots visible
  document.getElementById('boxplot-1').setAttribute('opacity', 1);
  document.getElementById('boxplot-2').setAttribute('opacity', 1);
  // Get data for bars
  percentiles = percentileMetrics[modelID][chainID][residueID];
  bAvg = percentiles[2];
  fitSC = percentiles[5];
  // Set main line coordinates
  bar1Y = parseFloat((barOffsetY + barMultiplierY * bAvg).toFixed(1));
  bar2Y = parseFloat((barOffsetY + barMultiplierY * fitSC).toFixed(1));
  document.getElementById('bar-1-mainline').setAttribute('y1', bar1Y);
  document.getElementById('bar-1-mainline').setAttribute('y2', bar1Y);
  document.getElementById('bar-2-mainline').setAttribute('y1', bar2Y);
  document.getElementById('bar-2-mainline').setAttribute('y2', bar2Y);
  // Set bar label text and position 
  let bar1Label = document.getElementById('bar-1-label');
  let bar2Label = document.getElementById('bar-2-label');
  bar1Label.textContent = bAvg;
  bar2Label.textContent = fitSC;
  if (bAvg < 10) {
    bar1Label.setAttribute('y', bar1Y-20);
  } else {
    bar1Label.setAttribute('y', bar1Y+20);
  };
  if (fitSC < 10) {
    bar2Label.setAttribute('y', bar2Y-20);
  } else {
    bar2Label.setAttribute('y', bar2Y+20);
  };
  // Set summary text
  let seqNum = sequenceNumbers[modelID][chainID][residueID];
  let aaCode = residueCodes[modelID][chainID][residueID];
  document.getElementById('residue-summary').textContent = 'Chain ' + String.fromCharCode(65 + chainID) + ', ' + 'Residue ' + seqNum + ' (' + aaCode + ')';
};


function clearResidueChart() {
  // Clear checkbox backgrounds
  document.getElementById('checkbox-1').setAttribute('fill', COLORS['VL_GREY']);
  document.getElementById('checkbox-2').setAttribute('fill', COLORS['VL_GREY']);
  // Clear checkbox text
  document.getElementById('checkbox-1-text').textContent = '';
  document.getElementById('checkbox-2-text').textContent = '';
  // Make box plots invisible
  document.getElementById('boxplot-1').setAttribute('opacity', 0);
  document.getElementById('boxplot-2').setAttribute('opacity', 0);
  // Clear summary text
  document.getElementById('residue-summary').textContent = 'Select a residue...';
};


function getResidueChartRanges() {
  for (var i=0; i<numModels; ++i) {
    let bAvgs = [ ];
    let fitSCs = [ ];
    for (var j=0; j<numChains; ++j) {
      for (var k=0; k<chainLengths[j]; ++k) {
        if (percentileMetrics[i][j][k] !== null) {
          let bAvg = percentileMetrics[i][j][k][2];
          if (bAvg !== null) {
            bAvgs.push(bAvg);
          };
          let fitSC = percentileMetrics[i][j][k][5];
          if (fitSC !== null) {
            fitSCs.push(fitSC);
          };
        };
      };
    };
    let bMean = mean(bAvgs);
    let bStd = standardDeviation(bAvgs);
    let bLow = Math.max(0, bMean-bStd);
    let bHigh = Math.min(100, bMean+bStd);
    let bMin = Math.min.apply(null, bAvgs);
    //let bMin = percentile(bAvgs, 0.05);
    //let bMin = Math.max(0, bMean-2*bStd);
    let bMax = Math.max.apply(null, bAvgs);
    //let bMax = percentile(bAvgs, 0.95);
    //let bMax = Math.min(100, bMean+2*bStd);
    let fitMean = mean(fitSCs);
    let fitStd = standardDeviation(fitSCs);
    let fitLow = Math.max(0, fitMean-fitStd);
    let fitHigh = Math.min(100, fitMean+fitStd);
    let fitMin = Math.min.apply(null, fitSCs);
    //let fitMin = percentile(fitSCs, 0.05);
    //let fitMin = Math.max(0, fitMean-2*fitStd);
    let fitMax = Math.max.apply(null, fitSCs);
    //let fitMax = percentile(fitSCs, 0.95);
    //let fitMax = Math.min(100, fitMean+2*fitStd);
    modelMinMax.push([ [bMean, bLow, bHigh, bMin, bMax], [fitMean, fitLow, fitHigh, fitMin, fitMax] ]);
  };
};


function setResidueChartRanges() {
  // Calculate Ys
  let bar1Mid = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][0][0]).toFixed(1));
  let bar1Low = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][0][1]).toFixed(1));
  let bar1High = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][0][2]).toFixed(1));
  let bar1Min = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][0][3]).toFixed(1));
  let bar1Max = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][0][4]).toFixed(1));
  let bar2Mid = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][1][0]).toFixed(1));
  let bar2Low = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][1][1]).toFixed(1));
  let bar2High = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][1][2]).toFixed(1));
  let bar2Min = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][1][3]).toFixed(1));
  let bar2Max = parseFloat((barOffsetY + barMultiplierY * modelMinMax[selectedModel][1][4]).toFixed(1));
  // Set box Ys
  let boxPlot1 = document.getElementById('boxplot-1-box');
  let boxPlot1Points = [ ];
  for (var i=0; i<boxPlot1.points.numberOfItems; ++i) {
    let x = boxPlot1.points.getItem(i).x;
    let y = boxPlot1.points.getItem(i).y;
    boxPlot1Points.push([x, y]);
  };
  boxPlot1Points[1][1] = boxPlot1Points[2][1] = bar1Min;
  boxPlot1Points[0][1] = boxPlot1Points[3][1] = bar1Max;
  boxPlot1.setAttribute('points', boxPlot1Points);
  let boxPlot2 = document.getElementById('boxplot-2-box');
  let boxPlot2Points = [ ];
  for (var i=0; i<boxPlot2.points.numberOfItems; ++i) {
    let x = boxPlot2.points.getItem(i).x;
    let y = boxPlot2.points.getItem(i).y;
    boxPlot2Points.push([x, y]);
  };
  boxPlot2Points[1][1] = boxPlot2Points[2][1] = bar2Min;
  boxPlot2Points[0][1] = boxPlot2Points[3][1] = bar2Max;
  boxPlot2.setAttribute('points', boxPlot2Points);
  // Set line Ys
  document.getElementById('boxplot-1-line-high').setAttribute('y1', bar1High);
  document.getElementById('boxplot-1-line-high').setAttribute('y2', bar1High);
  document.getElementById('boxplot-2-line-high').setAttribute('y1', bar2High);
  document.getElementById('boxplot-2-line-high').setAttribute('y2', bar2High);
  document.getElementById('boxplot-1-line-mid').setAttribute('y1', bar1Mid);
  document.getElementById('boxplot-1-line-mid').setAttribute('y2', bar1Mid);
  document.getElementById('boxplot-2-line-mid').setAttribute('y1', bar2Mid);
  document.getElementById('boxplot-2-line-mid').setAttribute('y2', bar2Mid);
  document.getElementById('boxplot-1-line-low').setAttribute('y1', bar1Low);
  document.getElementById('boxplot-1-line-low').setAttribute('y2', bar1Low);
  document.getElementById('boxplot-2-line-low').setAttribute('y1', bar2Low);
  document.getElementById('boxplot-2-line-low').setAttribute('y2', bar2Low);
  // Only show SD lines if they fall within the min/max range
  document.getElementById('boxplot-1-line-low').setAttribute('opacity', 0);
  if (bar1Low < bar1Min) {
    document.getElementById('boxplot-1-line-low').setAttribute('opacity', 1);
  }
  document.getElementById('boxplot-1-line-high').setAttribute('opacity', 0);
  if (bar1High > bar1Max) {
    document.getElementById('boxplot-1-line-high').setAttribute('opacity', 1);
  }
  document.getElementById('boxplot-2-line-low').setAttribute('opacity', 0);
  if (bar2Low < bar2Min) {
    document.getElementById('boxplot-2-line-low').setAttribute('opacity', 1);
  }
  document.getElementById('boxplot-2-line-high').setAttribute('opacity', 0);
  if (bar2High > bar2Max) {
    document.getElementById('boxplot-2-line-high').setAttribute('opacity', 1);
  }
  // Calculate and set gradient positions
  /*
  let gradient1PC = 1 - (modelMinMax[selectedModel][0][0] / (modelMinMax[selectedModel][0][4] - modelMinMax[selectedModel][0][3])).toFixed(3);
  let gradient2PC = 1 - (modelMinMax[selectedModel][1][0] / (modelMinMax[selectedModel][1][4] - modelMinMax[selectedModel][1][3])).toFixed(3);
  document.getElementById('gradient-1').childNodes[1].setAttribute('offset', gradient1PC);
  document.getElementById('gradient-2').childNodes[1].setAttribute('offset', gradient2PC);
  */
};


function main() {
  getResidueChartDims();
  getResidueChartRanges();
  setResidueChartRanges();
  if (numModels == 1) {
    document.getElementById('model-toggle').disabled = true;
  };
};


window.addEventListener('load', main);

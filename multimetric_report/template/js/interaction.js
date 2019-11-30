function scrollToId(elementId) {
	let el = document.getElementById(elementId);
	el.scrollIntoView({
		behavior: 'smooth',
		block: 'start'
	});
}

function toggleFullscreen() {
	if (chartZoomed) {
		chartZoomed = false;
		document.getElementById('radial-view').className = 'col-lg-6';
		document.getElementById('residue-view').className = 'col-lg-6';
		document.getElementById('radial-toggle-button').children[0].innerHTML = '↗';
	} else {
		chartZoomed = true;
		document.getElementById('radial-view').className = 'col-lg-8';
		document.getElementById('residue-view').className = 'col-lg-4';
		document.getElementById('radial-toggle-button').children[0].innerHTML = '↙';
	}
}

function setChain(chainID) {
	let radialCharts = document.querySelectorAll('[id^=radial-chart-]');
	for (var i = 0; i < radialCharts.length; ++i) {
		radialCharts[i].style.display = 'none';
		if (radialCharts[i].id === 'radial-chart-' + chainID) {
			radialCharts[i].style.display = '';
		}
	}
	let selectorButtons = document.querySelectorAll('[id^=chain-button-]');
	for (var i = 0; i < selectorButtons.length; ++i) {
		selectorButtons[i].style['color'] = null;
		if (selectorButtons[i].id === 'chain-button-' + chainID) {
			selectorButtons[i].style['color'] = 'rgb(150,150,150)';
		}
	}
}

function handleChain(actionID, chainID, residueID) {
	let residue = protein[chainID][residueID];
	let residueSeg = document.getElementById('radial-chart-' + chainID).getElementById('seg' + residueID);
	// On click
	if (actionID === 1) {
		document.getElementById('residue-summary').innerHTML = 'Chain ' + (chainID+1) + ', Residue ' + (residueID+1) + ' (' + residue.aa + ')';
		if (lastClicked) {
			let lastClickedSeg = document.getElementById('radial-chart-' + lastClicked.chainID).getElementById('seg' + lastClicked.residueID);
			lastClickedSeg.setAttribute('stroke-opacity', 0);
			lastClickedSeg.setAttribute('stroke-width', 1);
		}
		residueSeg.setAttribute('stroke-opacity', 1);
		residueSeg.setAttribute('stroke-width', 3);
		setRadarChart(residue.metrics);
		lastClicked = { 'chainID' : chainID, 'residueID' : residueID, 'residue' : residue }
	// On mouse over
	} else if (actionID === 2) {
		residueSeg.setAttribute('fill-opacity', 0.25);
	// On mouse out
	} else if (actionID === 3) {
		residueSeg.setAttribute('fill-opacity', 0);
		if (lastClicked) { document.getElementById('residue-summary').innerHTML = 'Chain ' + (lastClicked.chainID+1) + ', Residue ' + (lastClicked.residueID+1) + ' (' + lastClicked.residue.aa + ')'; }
	}
}

function handleSegment(actionID, chainID, residueID) {
	let residue = protein[chainID][residueID];
	let mainSegment = document.getElementById('segment' + chainID);
	mainSegment.setAttribute('transform', 'rotate(' + 360/protein[chainID].length*residueID + ', 500, 500' + ')');
}

function handleRadar(actionID, metricID) {
	let metricPoint = document.getElementById('metric-point-' + metricID);
	let floatingLabel = document.getElementById('floating-label-group');
	let floatingLabelText = document.getElementById('floating-label-text');
	// On mouse over
	if (actionID === 2) {
		floatingLabelText.innerHTML = lastClicked.residue.metrics[metricID];
		floatingLabel.setAttribute('transform', 'translate(' + (parseInt(metricPoint.getAttribute('cx'))) + ', ' + (parseInt(metricPoint.getAttribute('cy'))-30) + ')');
		floatingLabel.setAttribute('opacity', 1);
	// On mouse out
	} else if (actionID === 3) {
		floatingLabel.setAttribute('opacity', 0);
		floatingLabelText.innerHTML = '';
	}
}

function coordsFromAngle(centre, angle, pointRadius, adj=[0, 0]) {
	let xc = centre[0];
	let yc = centre[1];
	let x1 = xc + pointRadius * Math.sin(angle) + adj[0];
	let y1 = yc - pointRadius * Math.cos(angle) + adj[1];
	return [Math.round(x1, 2), Math.round(y1, 2)];
}

function setRadarChart(metricValues) {
	let radarChart = document.getElementById('radar-chart');
	let canvas = [parseInt(radarChart.getAttribute('viewBox').split(' ')[2]), parseInt(radarChart.getAttribute('viewBox').split(' ')[3])];
	let centre = [canvas[0]/2, canvas[1]/2];
	let angleDelta = 2*Math.PI / metricValues.length;
	let axisRadius = Math.min(canvas[0], canvas[1])/2 - Math.min(canvas[0], canvas[1])/20;
	let polygonPoints = [ ];
	for (var i = 0; i < metricValues.length; ++i) {
		let circle = document.getElementById('metric-point-' + i);
		let newCentre = coordsFromAngle(centre, i*angleDelta, axisRadius*metricValues[i]/100);
		polygonPoints.push(newCentre)
		circle.setAttribute('cx', newCentre[0]);
		circle.setAttribute('cy', newCentre[1]);
		let colour = metricValues[i] < 33 ? 'rgb(200, 50, 50)' : metricValues[i] < 67 ?'rgb(200, 200, 50)' : 'rgb(50, 200, 50)';
		circle.setAttribute('fill', colour);
		circle.setAttribute('fill-opacity', 0.5);
		circle.setAttribute('stroke-opacity', 1);
	};
	let polygon = document.getElementById('connecting-polygon');
	polygon.setAttribute('points', polygonPoints);
	polygon.setAttribute('fill-opacity', 0.5);
	polygon.setAttribute('stroke-opacity', 1);
}
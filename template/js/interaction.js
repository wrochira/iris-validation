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
		document.getElementById('chain-view').className = 'col-lg-7';
		document.getElementById('residue-view').className = 'col-lg-5';
		document.getElementById('toggle-button').children[0].innerHTML = '&#8599;';
	} else {
		chartZoomed = true;
		document.getElementById('chain-view').className = 'col-lg-8';
		document.getElementById('residue-view').className = 'col-lg-4';
		document.getElementById('toggle-button').children[0].innerHTML = '&#8601;';
	}
}

function setConcentric(chartID) {
	segmentDocked = true;
	let barCharts = document.querySelectorAll('[id^=bar-chart-]');
	for (var i = 0; i < barCharts.length; ++i) {
		barCharts[i].style.display = 'none';
		if (barCharts[i].id === 'bar-chart-' + chartID) {
			barCharts[i].style.display = '';
		}
	}
	let concentricCharts = document.querySelectorAll('[id^=concentric-chart-]');
	for (var i = 0; i < concentricCharts.length; ++i) {
		concentricCharts[i].style.display = 'none';
		if (concentricCharts[i].id === 'concentric-chart-' + chartID) {
			concentricCharts[i].style.display = '';
		}
	}
	let selectorButtons = document.querySelectorAll('[id^=chain-button-]');
	for (var i = 0; i < selectorButtons.length; ++i) {
		selectorButtons[i].style['color'] = null;
		if (selectorButtons[i].id === 'chain-button-' + chartID) {
			selectorButtons[i].style['color'] = 'rgb(150,150,150)';
		}
	}
}

function handleSegment(actionID, chainID, residueID, gapDegrees) {
	//let residue = protein[chainID][residueID];
	if (actionID === 1) {
		segmentDocked = false;
	}
	else if (actionID === 2 && !segmentDocked) {
		let mainSegment = document.getElementById('segment' + chainID);
		let startIndex = Math.floor(residueID-sidepanelHeight/2);
		let endIndex = Math.floor(residueID+sidepanelHeight-1-sidepanelHeight/2);
		if (startIndex >= 0 && endIndex < protein[chainID].length) {
			// Good.
		}
		else if (startIndex < 0) {
			startIndex = 0;
			endIndex = sidepanelHeight-1;
		}
		else {
			startIndex = protein[chainID].length-sidepanelHeight;
			endIndex = protein[chainID].length-1;
		}
		for (var i = 0; i < protein[chainID].length; i++) {
			let label = document.getElementById('residue-label-' + chainID + '-' + i);
			label.setAttribute('opacity', 0);
		};
		mainSegment.setAttribute('transform', 'rotate(' + (360-gapDegrees)/protein[chainID].length*startIndex + ', 500, 500' + ')');
		let startLabel = document.getElementById('residue-label-' + chainID + '-' + startIndex);
		let endLabel = document.getElementById('residue-label-' + chainID + '-' + endIndex.toString());
		startLabel.setAttribute('opacity', 1);
		endLabel.setAttribute('opacity', 1);
	}
	else if (actionID === 3) {
		segmentDocked = true;
	}
}

function handleRadar(actionID, metricID) {
	let metricPoint = document.getElementById('metric-point-' + metricID);
	let floatingLabel = document.getElementById('floating-label-group');
	let floatingLabelText = document.getElementById('floating-label-text');
	// On mouse over
	if (actionID === 2) {
		console.log(lastClicked)
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

function setRadarChart(chainID, residueID) {
	let residue = protein[chainID][residueID];
	let metricNames = Object.values(residue['metrics-absolute']);
	let metricValues = Object.values(residue['metrics-absolute']);
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
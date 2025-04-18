const {
	csv, 
	select, 
	scaleLinear, 
	scaleTime, 
	extent, 
	axisLeft, 
	axisBottom, 
	format,
	timeParse,
	line,
	max
} = d3;

// Set the dimensions and margins of the graph
const margin = {top: 10, right: 20, bottom: 50, left: 100};

function getContentWidth(element) {
	const styles = getComputedStyle(element);
	return (
		element.clientWidth - 
		parseFloat(styles.paddingLeft) -
		parseFloat(styles.paddingRight)
	);
}

function getDivWidth(element) {
	const styles = getComputedStyle(element);
	return (
		element.offsetWidth +
		parseFloat(styles.marginLeft) +
		parseFloat(styles.marginRight)
	);
}

function MBP(element) {
	const styles = getComputedStyle(element);
	return (
	parseFloat(styles.marginLeft) +
	parseFloat(styles.marginRight) +
	parseFloat(styles.borderLeft) +
	parseFloat(styles.borderRight) + 
	parseFloat(styles.paddingLeft) +
	parseFloat(styles.paddingRight)
	);
}
 
const width = 
	getContentWidth(document.getElementById('container')) - 
	getDivWidth(document.getElementById('leftChild')) -
	MBP(document.getElementById('leftChild')) -
	5;

const height = 
	window.innerHeight -
	MBP(document.getElementById('container')) -
	MBP(document.getElementById('leftChild')) -
	10;

// Append the svg object
const svg = select('#my_dataviz')
	.append('svg')
	.attr('width', width)
	.attr('height', height);

// Load and parse CSV
const csvUrl = "https://raw.githubusercontent.com/RyanPBlake/coal-project/main/coal-output-uk-tonnes.csv";
const parseRow = (d) => {
	d.Year = +d.Year;
	d.Coal_Output_BEIS_2020 = +d.Coal_Output_BEIS_2020;
	d.Coal_Imports_BEIS_2020 = +d.Coal_Imports_BEIS_2020;
	return d;
};

const xValue = (d) => d.Year;
const yValue = (d) => d.Coal_Output_BEIS_2020;
const xTime = (d) => d3.timeParse("%Y")(d.Year);
const radius = 3;
const commaFormat = format(',');

const coefficients = [29969978332, 28.48, 1921.74];

// Interpolate the data to fill missing years
function interpolate(data){
	const xyinterpol = [];
	const N = data.length;
	
	for (let i = 0; i < N - 1; i++) {
		xyinterpol.push({
			Year: data[i].Year,
			Coal_Output_BEIS_2020: data[i].Coal_Output_BEIS_2020
		});
		
		delta_X = data[i + 1].Year - data[i].Year;
		if (delta_X > 1) { 
			const delta_Y = 
				data[i + 1].Coal_Output_BEIS_2020 - data[i].Coal_Output_BEIS_2020;
			const m = delta_Y / delta_X;
			
			for (let year = 1; year < delta_X; year++) {
				const value = data[i].Coal_Output_BEIS_2020 + m * year;
				xyinterpol.push({
					Year : data[i].Year + year,
					Coal_Output_BEIS_2020 : value}
				);
			}
		} 
	}
	
	xyinterpol.push({
		Year : data[N-1].Year, 
		Coal_Output_BEIS_2020 : data[N-1].Coal_Output_BEIS_2020
	});
	
	return xyinterpol;
}

// Hubberts Curve Generator
function HubbertsCurve(coefficients, data) {
	const gamma = 1/coefficients[1];
	const points = [];
	
	for (let j = 0; j < data.length; j++) {
		const h_i = Math.E ** ((-gamma) * (data[j].Year - coefficients[2]));
		const g_i = coefficients[0] / (1 + h_i);
		const dgdx_i = g_i * (1 - (g_i / coefficients[0])) * gamma;
		points.push({Year : data[j].Year, estimate : dgdx_i});
	}
	
	return points;
}

// function that calculates the mean squared error between the estimate and the data
function MSE(estimate, data) {
	let residue = 0;
	const N = data.length;
	for (let j = 0; j < N; j++) {
		residue += (estimate[j].estimate - data[j].Coal_Output_BEIS_2020) ** 2;
	}
	return residue / N;
}

// === Gradient Descent Optimizer for Hubbert's Curve ===
async function optimizeHubbertCoefficients(data, initialCoefficients, learningRate = 1e-12, iterations = 1000) {
	const interpolatedData = interpolate(data);
	
	const loss = (coeffs) => {
		const curve = HubbertsCurve(coeffs, interpolatedData);
		return MSE(curve, interpolatedData);
	};
	
	function clamp(x, min, max) {
		return Math.min(Math.max(x, min), max);
	}
	
	const grad = (coeffs) => {
		const h = 10e-5;
		const grads = [0, 0, 0];
		
		// Central difference for each coefficient
		for (let i = 0; i < coeffs.length; i++) {
			const coeffsUp = [...coeffs];
			const coeffsDown = [...coeffs];
			
			// Adjust step size based on coefficient magnitude
            const step = h * Math.max(1, Math.abs(coeffs[i]));

			coeffsUp[i] += step;
			coeffsDown[i] -= step;
			
			const lossUp = loss(coeffsUp);
			const lossDown = loss(coeffsDown);
			
			grads[i] = (lossUp - lossDown) / (2 * step);
		}
		
		return grads;
	};
	
	let coeffs = [...initialCoefficients];
	let prevLoss = Infinity;
	
	for (let i = 0; i < iterations; i++) {
		const gradient = grad(coeffs);
		
		// Update coefficients with momentum
		coeffs = coeffs.map((c, idx) => {
			let updated = c - learningRate * gradient[idx];
			// Apply different constraints for each coefficient
			if (idx === 0) updated = clamp(updated, 1e7, 1e12);	// Q
			if (idx === 1) updated = clamp(updated, 1, 100);	// b
			if (idx === 2) updated = clamp(updated, 1900, 2100);	// t0
			return updated;
		});
		
		const currentLoss = loss(coeffs);
		const lossDiff = prevLoss - currentLoss;
		
		// Adaptive learning rate
		if (lossDiff > 0) {
			learningRate *= 1.05; // Increase if improving
        } else {
            learningRate *= 0.5;  // Decrease if not improving
        }
		
		prevLoss = currentLoss;
		
		if (i % 50 === 0) {
			console.log(`Iteration ${i}: MSE = ${loss(coeffs).toFixed(2)}`);
			console.log(`Coefficients: [${coeffs.map(c => c.toExponential(4)).join(', ')}]`);
			console.log(`Gradients: [${gradient.map(g => g.toExponential(4)).join(', ')}]`);
		}
		
		// Early stopping if improvement is minimal
		if (Math.abs(lossDiff) < 1e-6) {
			console.log(`Converged at iteration ${i}`);
			break;
		}
	}
	
	console.log("Optimized Coefficients:", coeffs);
	return coeffs;
}

const main = async function() {
	let data;
	try {
		data = await csv(csvUrl, parseRow);
	} catch (err) {
		console.error("Failed to load CSV:", err);
		return;
	}
	
	const interpolatedData = interpolate(data);
	
	const xRender = scaleTime()
		.domain(extent(data, xTime))
		.range([margin.left, width-margin.right]);
	
	const yRender = scaleLinear()
		.domain([0, max(data, yValue)])
		.range([height - margin.bottom, margin.top]);	
	
	// X axis	
	svg.append('g')
		.attr('transform', `translate(0,${height - margin.bottom})`)
		.call(axisBottom(xRender).tickFormat(d3.timeFormat('%Y')));
	
	// Y axis
	svg.append('g')
		.attr('transform', `translate(${margin.left},0)`)
		.call(axisLeft(yRender));

	const marks = data.map(d => ({
		x: xRender(xTime(d)), 
		y: yRender(yValue(d)),
		title : `Year: ${xValue(d)} | Output : ${commaFormat(yValue(d))}`
	}));

	// Line: Actual Data
    svg.append("path")
      	.datum(data)
      	.attr("fill", "none")
      	.attr("stroke", "steelblue")
      	.attr("stroke-width", 1.5)
      	.attr("d", line()
			.x(d => xRender(xTime(d)))
			.y(d => yRender(d.Coal_Output_BEIS_2020))
		);

	// Circles: Actual Data Points
	svg.selectAll('circle')
		.data(marks)
		.join('circle')
		.attr('cx', (d) => d.x)
		.attr('cy', (d) => d.y)
		.attr('r', radius)
		.each(function(d) {
			d3.select(this)
				.append('title')
				.text(d.title);
		});
	

	// X Label
	svg.append("text")
		.attr("text-anchor", "end")
		.attr("x", width / 2)
		.attr("y", height - 10)
		.text("Year");

	// Y Label
	svg.append("text")
		.attr("text-anchor", "end")
		.attr("transform", "rotate(-90)")
		.attr("x", -height / 2)
		.attr("y", 20)
		.text("Coal Output (Tonnes)");


	// Hubbert's Curve
	const optimizedCoefficients = await optimizeHubbertCoefficients(data, coefficients);
	const HubbertsPoints = HubbertsCurve(optimizedCoefficients, interpolatedData);
	svg.append("path")
		.datum(HubbertsPoints)
		.attr("fill", "none")
		.attr("stroke", "orange")
		.attr("stroke-width", 1.5)
		.attr("d", line()
			.x(d => xRender(xTime(d)))
			.y(d => yRender(d.estimate))
		);

	// MSE Display
	const MSEvalue =  MSE(HubbertsPoints, interpolatedData);
	select('#Text')	
		.append('p')
		.text(`MSE: ${commaFormat(MSEvalue)}`);
		
	// Legend
	const legendX = width - 150;
	const legendY = 30;

	svg.append("circle")
		.attr("cx", legendX)
		.attr("cy", legendY)
		.attr("r", 6)
		.style("fill", "steelblue");
	svg.append("text")
		.attr("x", legendX + 10)
		.attr("y", legendY)
		.text("Actual Output")
		.style("font-size", "12px")
		.attr("alignment-baseline", "middle");
		
	svg.append("circle")
		.attr("cx", legendX)
		.attr("cy", legendY + 20)
		.attr("r", 6).style("fill", "orange");
	svg.append("text")
		.attr("x", legendX + 10)
		.attr("y", legendY + 20)
		.text("Hubbert's Curve")
		.style("font-size", "12px")
		.attr("alignment-baseline", "middle");

};

main();




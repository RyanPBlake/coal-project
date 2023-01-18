const {csv, select, scaleLinear, scaleTime, extent, axisLeft, axisBottom, format} = d3;

// set the dimensions and margins of the graph
const margin = {top : 10, right: 20, bottom : 20, left : 100};

//function that return the content width of an element
function getContentWidth (element) {
	const styles = getComputedStyle(element)
  
	return element.clientWidth
	  - parseFloat(styles.paddingLeft)
	  - parseFloat(styles.paddingRight)
}

function getDivWidth (element) {
	const styles = getComputedStyle(element)

	return element.offsetWidth
		+ parseFloat(styles.marginLeft)
		+ parseFloat(styles.marginRight)
}

function MBP (element) {
	const styles = getComputedStyle(element)

	return parseFloat(styles.marginLeft)
		+ parseFloat(styles.marginRight)
		+ parseFloat(styles.borderLeft)
		+ parseFloat(styles.borderRight)
		+ parseFloat(styles.paddingLeft)
		+ parseFloat(styles.paddingRight)
}
 
const width = getContentWidth(document.getElementById('container')) - getDivWidth(document.getElementById('leftChild')) - MBP(document.getElementById('leftChild'))-5;
const height = window.innerHeight - MBP(document.getElementById('container')) - MBP(document.getElementById('leftChild'))-10;

// append the svg object to the specified div of the page
const svg = select('#my_dataviz')
	.append('svg')
	.attr('width', width)
	.attr('height', height);

// load the CSV Url and parse arguments into constants
const csvUrl = "https://raw.githubusercontent.com/RyanPBlake/coal-project/main/coal-output-uk-tonnes.csv";
const parseRow = function(d){
	d.Year = +d.Year;
	d.Coal_Output_BEIS_2020 = +d.Coal_Output_BEIS_2020;
	d.Coal_Imports_BEIS_2020 = +d.Coal_Imports_BEIS_2020;
	return d;
};
const xValue = function(d){return d.Year};
const yValue = function(d){return d.Coal_Output_BEIS_2020};
const xTime = function(d){return d3.timeParse("%Y")(d.Year)};
const radius = 3;
const commaFormat = format(',')


const coefficients = [29969978332, 28.48, 1921.74]

//interpolate the data to account for all the years
const interpolate = function(data){
	xyinterpol = [];
	N = data.length
	for (let i=0; i< N-1; i++){
		//push the current data
		xyinterpol.push({Year : data[i].Year, Coal_Output_BEIS_2020 : data[i].Coal_Output_BEIS_2020});
		
		//check to see if there is gaps in the data
		delta_X = data[i+1].Year - data[i].Year
		if (delta_X > 1) { 
			//create a linear interpolation y = mx+c that passes through (x_i,y_i) and (x_(i+1),y_(i+1))
			delta_Y = data[i+1].Coal_Output_BEIS_2020 - data[i].Coal_Output_BEIS_2020;
			m = delta_Y / delta_X;
			
			//so y = y_i + m*(x-x_i) but x is in years which is always an integer,  so x-x_i can be incremented as the variable "year"
			for (let year = 1; year < delta_X; year++){
				value = data[i].Coal_Output_BEIS_2020 + m*year
				xyinterpol.push({Year : data[i].Year + year, Coal_Output_BEIS_2020 : value});
			}
		} 
	}
	//push the final row of data
	xyinterpol.push({Year : data[N-1].Year, Coal_Output_BEIS_2020 : data[N-1].Coal_Output_BEIS_2020})
	return xyinterpol
}


// function that calculates Hubberts Curve from csv data and coefficients as input
const HubbertsCurve = function(coefficients, data){
	const gamma = 1/coefficients[1]
	points = []
	for (let j = 0; j < data.length; j++) {
		const h_i = Math.E ** ((-gamma) * (data[j].Year - coefficients[2]));
		const g_i = coefficients[0] / (1 + h_i);
		const dgdx_i = g_i * (1 - (g_i / coefficients[0])) * gamma;
		points.push({Year : data[j].Year, estimate : dgdx_i});
	}
	return points
};

// function that calculates the mean squared error between the estimate and the data
const MSE = function(estimate, data){
	let residue = 0;
	const N = data.length
	for (let j = 0; j < N; j++){
		residue += (estimate[j].estimate - data[j].Coal_Output_BEIS_2020)**2
	}
	residue = residue/N
	return residue
};



const main = async function(){
	// read the data and load it into a const
	const data = await csv(csvUrl, parseRow);
	const interpolatedData = interpolate(data)
	

	//Add x axis --> needs to be put into time format
	const xRender = scaleTime()
		.domain(extent(data, xTime))
		.range([margin.left, width-margin.right]);
	svg.append('g')
		.attr('transform',`translate(0,${height-margin.bottom})`)
		.call(axisBottom(xRender));


	//Add Y axis
	const yRender = scaleLinear()
		.domain([0, d3.max(data, yValue)])
		.range([height-margin.bottom, margin.top]);
	// y axis
	svg.append('g')
		.attr('transform',`translate(${margin.left},0)`)
		.call(axisLeft(yRender));

	//load pixel values of x and y coordinates into an object
	const marks = data.map(function(d){
		return {
			x : xRender(xTime(d)), 
			y : yRender(yValue(d)),
			title : `Year : ${xValue(d)} Output : ${commaFormat(yValue(d))})`,
		};
	});


	// Add the line on the raw data in blue
    svg
		.append("path")
      	.datum(data)
      	.attr("fill", "none")
      	.attr("stroke", "steelblue")
      	.attr("stroke-width", 1.5)
      	.attr("d", d3.line()
			.x(function(d){return xRender(xTime(d))})
			.y(function(d){return yRender(d.Coal_Output_BEIS_2020)})
		)

	//Add data points
	svg
		.selectAll('circle')
		.data(marks)
		.join('circle')
		.attr('cx', (d) => d.x)
		.attr('cy', (d) => d.y)
		.attr('r', radius)
		.append('title')
		.text(function(d){return d.title});

	


	//add hubbert's curve in orange
	const HubbertsPoints = HubbertsCurve(coefficients, interpolatedData);
	svg
		.append("path")
		.datum(HubbertsPoints)
		.attr("fill", "none")
		.attr("stroke", "orange")
		.attr("stroke-width", 1.5)
		.attr("d", d3.line()
			.x(function(d){return xRender(xTime(d))})
			.y(function(d){return yRender(d.estimate)})
		)

	//add the MES to the left hand side
	const MSEvalue =  MSE(HubbertsPoints, data);
	const MSErender = select('#Text')
						.append('p')
						.text(commaFormat(MSEvalue));

};

main();




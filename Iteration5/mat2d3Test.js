var margin = {top: 20, right: 80, bottom: 80, left: 80},
	width = 960 - margin.left - margin.right,
	height = 500 - margin.top - margin.bottom;

var xScale = d3.scale.linear()
	.domain([0, 7])
	.range([0, width]);

var yScale = d3.scale.linear()
	.domain([-1, 1])
	.range([height, 0]);

var xAxis = d3.svg.axis()
	.scale(xScale)
	.orient("bottom");

var yAxis = d3.svg.axis()
	.scale(yScale)
	.orient("left");

var line = d3.svg.line()
	.interpolate("linear")
	.x(function(d) { return xScale(d.x1); })
	.y(function(d) { return yScale(d.y1); });

var svg = d3.select("body").append("svg")
	.attr("width", width + margin.left + margin.right)
	.attr("height", height + margin.top + margin.bottom)
	.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

svg.append("g")
.attr("class", "axis")
.attr("transform", "translate(0, " + height + ")")
.call(xAxis);

svg.append("g")
.attr("class", "axis")
.call(yAxis);

d3.csv("mat2d3TestData1.csv", function(mydata)
{
	svg.append("path")
	.datum(mydata)
	.attr("class", "line1")
	.attr("d", line);

});

d3.csv("mat2d3TestData2.csv", function(mydata)
{
	svg.append("path")
	.datum(mydata)
	.attr("class", "line2")
	.attr("d", line);

});

d3.csv("mat2d3TestData3.csv", function(mydata)
{
	svg.append("path")
	.datum(mydata)
	.attr("class", "line3")
	.attr("d", line);

});

d3.csv("mat2d3TestData4.csv", function(mydata)
{
	svg.append("path")
	.datum(mydata)
	.attr("class", "line4")
	.attr("d", line);

});

d3.csv("mat2d3TestData5.csv", function(mydata)
{
	svg.append("path")
	.datum(mydata)
	.attr("class", "line5")
	.attr("d", line);

});


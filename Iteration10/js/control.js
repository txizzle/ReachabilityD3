var currTime = '1';
var currEvader = '1';
var contourDict = {};
var islands;
var leaveTrails = 0;
var showEvader = 1;
var evaderX;
var evaderY;

//Hash to map contours to number of islands
$.get("twoislands.txt", function(data, callback) {
  islands = data.split("\n");
});

//Initializing D3 necessities
var margin = {
      top: 80,
      right: 80,
      bottom: 80,
      left: 80
    },
    width = 500 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;
var xScale = d3.scale.linear().domain([1, 31]).range([0, width]);
var yScale = d3.scale.linear().domain([1, 31]).range([height, 0]);
var exScale = d3.scale.linear().domain([-1, 1]).range([0, width]);
var eyScale = d3.scale.linear().domain([-1, 1]).range([height, 0]);

var xAxis = d3.svg.axis().scale(xScale).orient("bottom");
var yAxis = d3.svg.axis().scale(yScale).orient("left");
var line = d3.svg.line().interpolate("linear").x(function(d) {
    return xScale(d.y1);
  }).y(function(d) {
    return yScale(d.x1);
  });

var objectLine = d3.svg.line().interpolate("linear").x(function(d) {
  return exScale(d[0]);
}).y(function(d) {
  return eyScale(d[1]);
});

//Initialize SVG
var svg = d3.select("body").append("svg").attr("width", width + margin.left +
  margin.right).attr("height", height + margin.top + margin.bottom).append(
  "g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
svg.append("g").attr("class", "axis").attr("transform", "translate(0, " +
  height + ")").call(xAxis);
svg.append("g").attr("class", "axis").call(yAxis);

$(document).ready(function() {
  //Plot contour for t=1, z=1
  d3.csv("../Iteration6/csv/t1z1Data1.csv", function(mydata) {
    svg.append("path").datum(mydata).attr("class", "line").attr("d", line);
  });

  d3.csv("../Iteration6/csv/t1z1Data2.csv", function(mydata) {
    svg.append("path").datum(mydata).attr("class", "line").attr("d", line);
  });

  //Plot evader path
  evaderX = getEvaderX(1);
  evaderY = getEvaderY(1);
  svg.append("path").datum([
    [-0.6, -0.8],
    [-0.6, 0.8],
    [0.6, 0.8],
    [0.6, -0.8],
    [-0.6, -0.8]
  ]).attr("class", "evaderpath").attr("d", objectLine).attr("stroke-dasharray",
    "10 5");

  //Plot Evader
  var evader = svg.append("circle").datum([evaderX, evaderY]).attr("cx", function(
    d) {
    return exScale(d[0]);
  }).attr("cy", function(d) {
    return eyScale(d[1]);
  }).attr("r", 5).attr("fill", "blue").call(drag);

  //Plot catch radius
  var catch_radius = svg.append("circle").datum([evaderX, evaderY]).attr("cx",
    function(d) {
      return exScale(d[0]);
    }).attr("cy", function(d) {
    return eyScale(d[1]);
  }).attr("r", height / 20).attr("fill", "none").attr("stroke", "black").attr(
    "stroke-width", "3px").attr("stroke-dasharray", "10 5");
});

function updateTime(t) {
  currTime = t.toString();
  $('#timeLabel').val(t);
  document.getElementById('timeSlider').value = t;
  updateContour();
}

function updateEvader(v) {
  currEvader = v.toString();
  $('#evaderLabel').val(v);
  evaderX = getEvaderX(v);
  evaderY = getEvaderY(v);
  document.getElementById('evaderSlider').value = v;
  updateContour();
}

$("input[name=optradio]:radio").change(function() {
  leaveTrails = $(this).val();
  d3.selectAll(".line").remove();
});

element.addEventListener("mousedown", function(e) {
  e.preventDefault();
}, false);
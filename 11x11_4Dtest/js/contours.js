var currTime = '1';
var currEvader = '1';
var currDefender = '1';
var contourDict = {};
var islands;
var leaveTrails = 0;
var special;

//Map from contours to number of 'islands' in plot
$.get("../data_contours/twoislands.txt", function(data) {
  islands = data.split("\n");
});

//Initial D3 Necessities
var margin = {
    top: 40,
    right: 40,
    bottom: 40,
    left: 40
  },
  width = 550 - margin.left - margin.right,
  height = 550 - margin.top - margin.bottom;
var xScale = d3.scale.linear().domain([1, 11]).range([0, width]);
var yScale = d3.scale.linear().domain([1, 11]).range([height, 0]);
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
var svg = d3.select("#contour").append("svg").attr("width", width + margin.left +
  margin.right).attr("height", height + margin.top + margin.bottom).append(
  "g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
svg.append("g").attr("class", "axis").attr("transform", "translate(0, " +
  height + ")").call(xAxis);
svg.append("g").attr("class", "axis").call(yAxis);

//Initial Contour for t=1, d=1, z=1
d3.csv("../data_contours/t1d1z2Data.csv", function(mydata) {
  svg.append("path").datum(mydata).attr("class", "line").attr("d", line);
});

function updateContour() {
  var currFile = "../data_contours/t" + currTime + "d" + currDefender + "z" + currEvader + "Data";
  var key = "t" + currTime + +"d" + currDefender + "z" + currEvader;
  d3.selectAll("path").attr("class", "line");
  if (leaveTrails == 0) {
    d3.selectAll("path").remove();
  }
  if (islands.indexOf(key) > -1) {
    d3.csv(currFile + "1.csv", function(mydata) {
      svg.append("path").datum(mydata).attr("class", "mainline").attr("d",
        line);
    });
    d3.csv(currFile + "2.csv", function(mydata) {
      svg.append("path").datum(mydata).attr("class", "mainline").attr("d",
        line);
    });
  } else {
    d3.csv(currFile + ".csv", function(mydata) {
      svg.append("path").datum(mydata).attr("class", "mainline").attr("d",
        line);
    });
  }
};

function updateTime(t) {
  currTime = t.toString();
  $('#timeLabel').val(t);
  document.getElementById('timeSlider').value = t;
  updateContour();
}
function updateEvader(v) {
  currEvader = v.toString();
  $('#evaderLabel').val(v);
  document.getElementById('evaderSlider').value = v;
  updateContour();
}
function updateDefender(v) {
  currDefender = v.toString();
  $('#defenderLabel').val(v);
  document.getElementById('defenderSlider').value = v;
  updateContour();
}


$("input[name=optradio]:radio").change(function() {
  leaveTrails = $(this).val();
  d3.selectAll(".line").remove();
});

element.addEventListener("mousedown", function(e) {
  e.preventDefault();
}, false);
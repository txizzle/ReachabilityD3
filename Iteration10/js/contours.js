var currTime = '1';
var currEvader = '1';
var contourDict = {};
var islands;
var leaveTrails = 0;
var evaderX;
var evaderY;

//Map from contours to number of 'islands' in plot
$.get("twoislands.txt", function(data) {
  islands = data.split("\n");
});

//Initial D3 Necessities
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

//Drag behavior for Evader
var drag = d3.behavior.drag().on('dragstart', function() {
  evader.style('fill', 'red');
}).on('drag', function() {
  max_x = exScale(0.6); //right
  min_x = exScale(-0.6); //left
  max_y = eyScale(-0.8); //bottom
  min_y = eyScale(0.8); //top
  new_x = exScale.invert(d3.event.x);
  new_y = eyScale.invert(d3.event.y);
  if (evaderX == 0.6) {
    //we are on the right side
    console.log("we on the right");
    if (evaderY == 0.8 || evaderY == -0.8) {
      //we are on a right corner
      console.log("we on a right corner");
      debugger;
      if (new_x < evaderX) { //moving left
        evaderX = Math.max(new_x, -0.6);
      } else {
        if (new_y > evaderY) { //moving up
          evaderY = Math.min(new_y, 0.8);
        } else if (new_y < evaderY) {
          console.log("moving down"); //moving down
          evaderY = Math.max(new_y, -0.8);
        }
      }
    } else {
      if (new_y > evaderY) { //moving up
        evaderY = Math.min(new_y, 0.8);
      } else if (new_y < evaderY) { // moving down
        evaderY = Math.max(new_y, -0.8);
      }
    }
  } else if (evaderX == -0.6) {
    //we are on the left side
    console.log("we on the left");
    if (evaderY == 0.8 || evaderY == -0.8) {
      //we are on a right corner
      console.log("we on a left corner");
      if (new_x > evaderX) { //moving right
        evaderX = Math.min(new_x, 0.6);
      } else {
        if (new_y > evaderY) { //moving up
          evaderY = Math.min(new_y, 0.8);
        } else if (new_y < evaderY) {
          console.log("moving down"); //moving down
          evaderY = Math.max(new_y, -0.8);
        }
      }
    } else {
      if (new_y > evaderY) { //moving up
        evaderY = Math.min(new_y, 0.8);
      } else if (new_y < evaderY) {
        console.log("moving down"); //moving down
        evaderY = Math.max(new_y, -0.8);
      }
    }
  } else if (evaderY == 0.8) {
    //we are on the top side
    console.log("we on the top");
    if (new_x > evaderX) { //moving right
      evaderX = Math.min(new_x, 0.6);
    } else if (new_x < evaderX) { //moving left
      evaderX = Math.max(new_x, -0.6);
    }
  } else if (evaderY == -0.8) {
    //we are on the bottom side
    console.log("we on the bottom");
    if (new_x > evaderX) { //moving right
      evaderX = Math.min(new_x, 0.6);
    } else if (new_x < evaderX) { //moving left
      evaderX = Math.max(new_x, -0.6);
    }
  }
  evader.attr('cx', exScale(evaderX)).attr('cy', eyScale(evaderY));
  catch_radius.attr('cx', exScale(evaderX)).attr('cy', eyScale(evaderY));
  currEvader = getValfromXY(evaderX, evaderY);
  //if we want to update contour here:
  //                 updateEvader(currEvader);
}).on('dragend', function() {
  evader.style('fill', 'black');
  updateEvader(currEvader);
});

//Initialize SVG
var svg = d3.select("body").append("svg").attr("width", width + margin.left +
  margin.right).attr("height", height + margin.top + margin.bottom).append(
  "g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
svg.append("g").attr("class", "axis").attr("transform", "translate(0, " +
  height + ")").call(xAxis);
svg.append("g").attr("class", "axis").call(yAxis);

//Initial Contour for t=1, z=1
d3.csv("../Iteration6/csv/t1z1Data1.csv", function(mydata) {
  svg.append("path").datum(mydata).attr("class", "line").attr("d", line);
});

d3.csv("../Iteration6/csv/t1z1Data2.csv", function(mydata) {
  svg.append("path").datum(mydata).attr("class", "line").attr("d", line);
});

evaderX = getEvaderX(1);
evaderY = getEvaderY(1);

//Plot evader path
svg.append("path").datum([
  [-0.6, -0.8],
  [-0.6, 0.8],
  [0.6, 0.8],
  [0.6, -0.8],
  [-0.6, -0.8]
]).attr("class", "evaderpath").attr("d", objectLine).attr("stroke-dasharray",
  "10 5");

//Plot evader
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

function updateContour() {
  var currFile = "../Iteration6/csv/t" + currTime + "z" + currEvader + "Data";
  var key = "t" + currTime + "z" + currEvader;
  d3.selectAll("path").attr("class", "line");
  if (leaveTrails == 0) {
    d3.selectAll("path").remove();
  }
    d3.selectAll("circle").remove();
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
  //Plot evader path
  svg.append("path").datum([
    [-0.6, -0.8],
    [-0.6, 0.8],
    [0.6, 0.8],
    [0.6, -0.8],
    [-0.6, -0.8]
  ]).attr("class", "evaderpath").attr("d", objectLine).attr(
    "stroke-dasharray", "10 5");
  
  //Plot evader
  evader = svg.append("circle").datum([evaderX, evaderY]).attr("cx", function(
    d) {
    return exScale(d[0]);
  }).attr("cy", function(d) {
    return eyScale(d[1]);
  }).attr("r", 5).attr("fill", "blue").call(drag);
  
  //Plot catch radius
  catch_radius = svg.append("circle").datum([evaderX, evaderY]).attr("cx",
    function(d) {
      return exScale(d[0]);
    }).attr("cy", function(d) {
    return eyScale(d[1]);
  }).attr("r", height / 20).attr("fill", "none").attr("stroke", "black").attr(
    "stroke-width", "3px").attr("stroke-dasharray", "10 5").call(drag);
  
  //Plot obstacle
 svg.append("path").datum(getObstacle(parseInt(currTime))).attr("class",
    "obstacle").attr("d", objectLine);
};

function updateTime(t) {
  currTime = t.toString();
  $('#timeLabel').val(t);
  document.getElementById('timeSlider').value = t;
  updateContour();
}

function getValfromXY(x, y) {
  var s;
  if (x == 0.6) {
    // 1.2 < s < 2.8
    s = y + 2;
  } else if (x == -0.6) {
    // 4 < s < 5.6
    s = 4.8 - y;
  } else {
    if (y == 0.8) { //top
      // 1.2 < s < 2.8
      s = 3.4 - x;
    } else if (y == -0.8) { //bottom
      // 4 < s < 5.6
      s = x + 0.6;
    }
  }
  return Math.round(s * 85 / 5.6 + 1);
}

function getEvaderY(val) {
  var s = (val - 1) / 85 * 5.6;
  if (s < 1.2) return -0.8;
  else if (s < 2.8) return -0.8 + (s - 1.2);
  else if (s < 4.0) return 0.8;
  else return 0.8 - (s - 4.0);
}

function getEvaderX(val) {
  var s = (val - 1) / 85 * 5.6;
  if (s < 1.2) return s - 0.6;
  else if (s < 2.8) return 0.6;
  else if (s < 4.0) return 0.6 - (s - 2.8);
  else return -0.6;
}

function getObstacle(time) {
  var t = time / 127.0;
  if (t >= 0 && t < 0.2) {
    var bt = -2 * t;
  } else if (t < 0.6) {
    var bt = -0.4;
  } else if (t < 0.8) {
    var bt = -2.8 + 4 * t;
  } else {
    var bt = 0.4;
  }
  return ([
    [-0.2, -0.6 + bt],
    [-0.2, 0.6 + bt],
    [0.2, 0.6 + bt],
    [0.2, -0.6 + bt]
  ]);
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
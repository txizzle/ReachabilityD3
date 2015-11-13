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

function updateContour() {
  var currFile = "../Iteration6/csv/t" + currTime + "z" + currEvader + "Data";
  var key = "t" + currTime + "z" + currEvader;
  d3.selectAll("path").attr("class", "line");
  
  if (leaveTrails == 0) {
    d3.selectAll("path").remove();
  }
  
  if (showEvader == 1) {
    d3.selectAll("circle").remove();
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
  
  //evader path
  svg.append("path").datum([
    [-0.6, -0.8],
    [-0.6, 0.8],
    [0.6, 0.8],
    [0.6, -0.8],
    [-0.6, -0.8]
  ]).attr("class", "evaderpath").attr("d", objectLine).attr(
    "stroke-dasharray", "10 5");
  
  //evader
  evader = svg.append("circle").datum([evaderX, evaderY]).attr("cx", function(
    d) {
    return exScale(d[0]);
  }).attr("cy", function(d) {
    return eyScale(d[1]);
  }).attr("r", 5).attr("fill", "blue").call(drag);
  //catch radius
  catch_radius = svg.append("circle").datum([evaderX, evaderY]).attr("cx",
    function(d) {
      return exScale(d[0]);
    }).attr("cy", function(d) {
    return eyScale(d[1]);
  }).attr("r", height / 20).attr("fill", "none").attr("stroke", "black").attr(
    "stroke-width", "3px").attr("stroke-dasharray", "10 5").call(drag);
  
  //obstacle
  svg.append("path").datum(getObstacle(parseInt(currTime))).attr("class",
    "obstacle").attr("d", objectLine);
  if (leaveTrails == 1) {
    d3.select(".line").remove();
  }
};

function updateTime(t) {
  currTime = t.toString();
  $('#timeLabel').val(t);
  document.getElementById('timeSlider').value = t;
  updateContour();
}

function getValfromXY(x, y) {
  var s;
  if (x == 0.6) { //right
    s = y + 2;
  } else if (x == -0.6) { //left
    s = 4.8 - y;
  } else {
    if (y == 0.8) { //top
      s = 3.4 - x;
    } else if (y == -0.8) { //bottom
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

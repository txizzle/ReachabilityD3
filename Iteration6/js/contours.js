var currTime = '2';
var currEvader = '1';
var currLevel = 0;
var contourDict = {};
var islands;
$.get("twoislands.txt", 
    function(data) {
        islands = data.split("\n");
        debugger;
    }
);

//initial svg contour
var margin = {top: 20, right: 80, bottom: 80, left: 80},
width = 500 - margin.left - margin.right,
height = 500 - margin.top - margin.bottom;

var xScale = d3.scale.linear()
    .domain([0, 30])
    .range([0, width]);

var yScale = d3.scale.linear()
    .domain([0, 30])
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

d3.csv("csv/t1z1Data1.csv", function(mydata)
{
    svg.append("path")
    .datum(mydata)
    .attr("class", "line")
    .attr("d", line);

});

d3.csv("csv/t1z1Data2.csv", function(mydata)
{
    svg.append("path")
    .datum(mydata)
    .attr("class", "line")
    .attr("d", line);

});
    //
    
//    var data = d3.csv.parseRows(text).map(function(row) {
//        return row.map(function(value) {
//            return -value;
//        });
//    });
//    var test = d3.select("#contour").append("svg")
//        .attr("width", width)
//        .attr("height", height)
//        .selectAll("path")
//        .data(c.contourList())
//        .enter().append("path")
//        .style("fill", function(d) {
//            return colours(d.level);
//        })
//        .style("stroke", "black")
//        .attr("d", d3.svg.line()
//            .x(function(d) {
//                return x(d.x);
//            })
//            .y(function(d) {
//                return y(d.y);
//            }));

function updateContour(i) {
    var key = currTime.toString() + "," + currEvader.toString() + "Data";
    //if (key in contourDict) {
    if (1 == 2) { //temporariliy disable memoization
        var c = contourDict[key];
        console.log("blahblah");
        //console.log(c);
        d3.select("svg").remove();
        var test = d3.select("#contour").append("svg")
            .attr("width", 300)
            .attr("height", 300)
            .selectAll("path")
            .data(c.contourList())
            .enter().append("path")
            //.style("fill",function(d) { return colours(d.level);})
            .style("stroke", "black")
            .attr("d", d3.svg.line()
                .x(function(d) {
                    return x(d.x);
                })
                .y(function(d) {
                    return y(d.y);
                }));
    } else {
        currLevel = i;
        var currFile = "csv/t" + currTime + "z" + currEvader +"Data";
        key = "t" + currTime + "z" + currEvader;
        d3.select("svg").remove();
        
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
        
        if (islands.indexOf(key) > -1) {
            console.log(">1 islands");
            d3.csv(currFile+"1.csv", function(mydata)
                {
                    svg.append("path")
                    .datum(mydata)
                    .attr("class", "line")
                    .attr("d", line);

                });

            d3.csv(currFile+"2.csv", function(mydata)
            {
                svg.append("path")
                .datum(mydata)
                .attr("class", "line")
                .attr("d", line);

            });
          }
          else{
            console.log("only one island");
            d3.csv(currFile+".csv", function(mydata)
            {
                svg.append("path")
                .datum(mydata)
                .attr("class", "line")
                .attr("d", line);

            });
            
          }
        };
};

var i = 0;

function updateTime(t) {
    currTime = t.toString();
    $('#timeLabel').val(t);
    document.getElementById('timeSlider').value = t;
    updateContour(currLevel);
}

function updateEvader(v) {
    currEvader = v.toString();
    $('#evaderLabel').val(v);
    document.getElementById('evaderSlider').value = v;
    updateContour(currLevel);
}

function loadJSON() {
    var json;
    $.getJSON("../value_json2.txt", function(data) {
        var items = [];
        $.each(data, function(key, val) {
            console.log(key);
            console.log(val);
            items.push(val);
            localStorage.setItem(key, val._ArrayData);
            debugger;
            $.cookie("test", val);
        });
    });
}
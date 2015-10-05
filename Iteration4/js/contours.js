var currTime = '2';
var currEvader = '1';
var currLevel = 0;
var contourDict = {};

//initial svg contour
d3.text("csv/t2z1.csv", function(text) {
    var data = d3.csv.parseRows(text).map(function(row) {
        return row.map(function(value) {
            return -value;
        });
    });
    //console.log(data);
    var cliff = -1000;
    data.push(d3.range(data[0].length).map(function() {
        return cliff;
    }));
    data.unshift(d3.range(data[0].length).map(function() {
        return cliff;
    }));
    data.forEach(function(d) {
        d.push(cliff);
        d.unshift(cliff);
    });

    var c = new Conrec,
        xs = d3.range(0, data.length),
        ys = d3.range(0, data[0].length),
        zs = [currLevel],
        width = 300,
        height = 300,
        x = d3.scale.linear().range([0, width]).domain([0, data.length]),
        y = d3.scale.linear().range([height, 0]).domain([0, data[0].length]),
        colours = d3.scale.linear().domain([-5, 3]).range(["#fff",
            "red"
        ]);

    c.contour(data, 0, xs.length - 1, 0, ys.length - 1, xs, ys, zs.length,
        zs);
    var test = d3.select("#contour").append("svg")
        .attr("width", width)
        .attr("height", height)
        .selectAll("path")
        .data(c.contourList())
        .enter().append("path")
        .style("fill", function(d) {
            return colours(d.level);
        })
        .style("stroke", "black")
        .attr("d", d3.svg.line()
            .x(function(d) {
                return x(d.x);
            })
            .y(function(d) {
                return y(d.y);
            }));
});

function updateContour(i) {
    var key = currTime.toString() + "," + currEvader.toString();
    //if (key in contourDict) {
    if (1 == 2) { //temporariliy disable memoization
        var c = contourDict[key];
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
        var currFile = "csv/t" + currTime + "z" + currEvader + ".csv";
        d3.text(currFile, function(text) {
            var data = d3.csv.parseRows(text).map(function(row) {
                return row.map(function(value) {
                    return -value;
                });
            });
            // console.log(data);
            var cliff = -1000;
            data.push(d3.range(data[0].length).map(function() {
                return cliff;
            }));
            data.unshift(d3.range(data[0].length).map(function() {
                return cliff;
            }));
            data.forEach(function(d) {
                d.push(cliff);
                d.unshift(cliff);
            });

            var c = new Conrec,
                xs = d3.range(0, data.length),
                ys = d3.range(0, data[0].length),
                zs = [i],
                width = 300,
                height = 300,
                x = d3.scale.linear().range([0, width]).domain([0,
                    data.length
                ]),
                y = d3.scale.linear().range([height, 0]).domain([0,
                    data[0].length
                ]);
            
            c.contour(data, 0, xs.length - 1, 0, ys.length - 1, xs,
                ys, zs.length, zs);
            
            contourDict[key] = c;
            d3.select("svg").remove();
            var test = d3.select("#contour").append("svg")
                .attr("width", width)
                .attr("height", height)
                .selectAll("path")
                .data(c.contourList())
                .enter().append("path")
                .style("stroke", "black")
                .attr("d", d3.svg.line()
                    .x(function(d) {
                        return x(d.x);
                    })
                    .y(function(d) {
                        return y(d.y);
                    }));
        });
    };
}

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
$(document).ready(function () {
    $("#roundslider").roundSlider({
        min: 1,
        max: 86,
        value: 1,
        radius: 50,
        drag: function (e) {
            updateEvader(e.value);
        },
        change: function (e) {
            updateEvader(e.value);
        }
    });
});
            
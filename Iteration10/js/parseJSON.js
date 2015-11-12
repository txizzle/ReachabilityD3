function loadJSON() {
    
$.getJSON( "../value_json2.txt", function( data ) {
  var items = [];
  $.each( data, function( key, val ) {
      console.log(key);
    console.log(val);
  });
});
}

var bedpe_feature = function(board){
// bedpe feature
  var bedpe_feature = tnt.board.track.feature();
    
  // Create
  bedpe_feature.create (function (elems) {
    var track = this;
    var xScale = bedpe_feature.scale();

    var y = track.height();
    //var x = track.width();
    var ymax = (board.width() / (xScale(10000) - xScale(0)) * 10000);
    var yScale = function(ycoord){
      return (y - (0 + ycoord)/ymax * y);
    };
    
    var format_xy = function(a, b){
      //console.log(xScale((a + b) / 2) + "," + (y - (2 + Math.abs(a - b))/4 * y));
      return xScale((a + b) / 2) + "," + yScale(Math.abs(a - b));
    };
 
    var d2points = function(d){
      if(d.s1 == d.s2 && d.e1 == d.e2){ //catch self interaction
        out = format_xy(d.s1, d.s2) + " " +
              format_xy(d.e1, d.e2) + " " +
              format_xy(d.e1, d.s2);
      }else{
        out = format_xy(d.s1, d.s2) + " " +
              format_xy(d.s1, d.e2) + " " +
              format_xy(d.e1, d.e2) + " " +
              format_xy(d.e1, d.s2);
      }
      return out;
    };
 
    var g = elems
    .append("g");
    g
    .append("polygon")
    .attr("points", function(d){
      return d2points(d);
    })
    .attr("fill", function(d){
      return d.color;
    });

  });
  
  // Move
  bedpe_feature.move (function (blocks) {
    var xScale = bedpe_feature.scale();
    var track = this;
    var y = track.height();
    //syncs yScale with xScale
    var ymax = (board.width() / (xScale(10000) - xScale(0)) * 10000);//board.to() - board.from();
    var yScale = function(ycoord){
      return (y - (0 + ycoord)/ymax * y);
    };
    var format_xy = function(a, b){
      return xScale((a + b) / 2) + "," + yScale(Math.abs(a - b));
    };
 
    var d2points = function(d){
      if(d.s1 == d.s2 && d.e1 == d.e2){ //catch self interaction
        out = format_xy(d.s1, d.s2) + " " +
              format_xy(d.e1, d.e2) + " " +
              format_xy(d.e1, d.s2);
      }else{
        out = format_xy(d.s1, d.s2) + " " +
              format_xy(d.s1, d.e2) + " " +
              format_xy(d.e1, d.e2) + " " +
              format_xy(d.e1, d.s2);
      }
      return out;
    };

   blocks.select("polygon")
    .attr("points", function(d){
      return d2points(d);
    });
    
  });
  
  return bedpe_feature;
};
  
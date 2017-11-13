var custom_feature = function (div) {
  
  var board = tnt.board().from(-1).to(5).min(-8).max(10).width(500).zoom_out(10).zoom_in(1);



  // Axis track
  var axis_track = tnt.board.track()
  .height(0)
  .color("white")
  .display(tnt.board.track.feature.axis()
           .orientation("top")
  );
  
  /*
  board.pane
  .append("text")
  .text("debug");
  */
  
  // arrow feature
  var polygon_feature = tnt.board.track.feature();
  
  // Create
  polygon_feature.create (function (elems) {
    
   
    var track = this;
    var xScale = polygon_feature.scale();
    var xresize = function(size){
      return board.width() * (size / (board.to() - board.from()));
    };
    var y = track.height();
    //var x = track.width();
 
    var g = elems
    .append("g")
    .attr("transform", function (d) {
      return "translate(" + xScale(d.mid) + "," + (y - (0.2 + d.dist)/(1 + 2)*y) + ")";
      //return "translate(" + (0 + d.mid) + "," + (y - (0.2 + d.dist)/(1 + 2)*y) + ")";
    });
    g
    .append("polygon")
    .attr("points", function(d){
      return -xresize(d.size) + ",0 0,25 " + xresize(d.size) + ",0 0,-25";
      //return "-25,0 0,25 25,0 0,-25";
    })
    .attr("fill", function(d){
      return d.color;
    });
    /*
    g
    .append("text")
    .attr("x", 0)
    .attr("y", 0)
    .attr("style", "fill:blue")
    .text("test");
    */

  });
  
  
  // Move
  polygon_feature.move (function (blocks) {
    var xScale = polygon_feature.scale();
    var track = this;
    var y = track.height();

/*
    blocks.select("polygon")
    .attr("points", function(d){
      return -xresize(d.size) + ",0 0,25 " + xresize(d.size) + ",0 0,-25";
      //return "-25,0 0,25 25,0 0,-25";
    });
  */  
       blocks.select("polygon")
    .attr("points", function(d){
      return xScale(board.from()) + ",0 0,25 " + xScale(board.to()) + ",0 0,-25";
      //return -board.width()/2 + ",0 0,25 " + board.width()/2 + ",0 0,-25";
      //return "-25,0 0,25 25,0 0,-25";
    });
    
  });
  

  
  
  // Data track
  var polygon_track = tnt.board.track()
  .height(60)
  .color("white")
  .display(polygon_feature)
  .data(tnt.board.track.data.sync()
        .retriever (function () {
          return [
            {
              mid : .5,
              size : 1,
              dist : 0,
              color : "red",
          
            },
            {
              mid : 2,
              size : 1,
              dist : 0,
              color : "orange",
          
            },
            {
              mid : 3,
              size : 1,
              dist : 0,
              color : "blue",
          
            },
            {
              mid : 1.5,
              size : 1,
              dist : 1,
              color : "darkred",
          
            },
            {
              mid : 2.5,
              size : 1,
              dist : 1,
              color : "darkorange",
          
            },
            {
              mid : 2,
              size : 1,
              dist : 2,
              color : "green",
          
            }
            ];
        })
  );
  
  
  board
  .add_track([axis_track, polygon_track]);
  
  board(div);
  
  board.start();
};

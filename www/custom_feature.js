var custom_feature = function (div) {
  
  var board = tnt.board().from(-1).to(5).min(-8).max(10).width(500).zoom_out(10).zoom_in(1);
  
  // Axis track
  var axis_track = tnt.board.track()
  .height(0)
  .color("white")
  .display(tnt.board.track.feature.axis()
           .orientation("top")
  );
  
  // arrow feature
  var arrow_feature = tnt.board.track.feature();
  
  // Create
  arrow_feature.create (function (elems) {
    var xScale = arrow_feature.scale();
   
    var track = this;
    var y = track.height();
    //var x = track.width();
 
    var g = elems
    .append("g")
    .attr("transform", function (d) {
      return "translate(" + xScale(d.mid) + "," + (y - (0.2 + d.dist)/(1 + 2)*y) + ")";
      //return "translate(" + (0 + d.mid) + "," + (y - (0.2 + d.dist)/(1 + 2)*y) + ")";
    });
    g
    .append("line")
    .attr("x1", 0)
    .attr("y1", 0)
    .attr("x2", 0)
    .attr("y2", function(d){
      return 0;
      //return -(1 + d.dist)/(1 + 2)*y;
    })
    .attr("stroke", function (d) {
      return d.color;
    });
    
    g
    .append("line")
    .attr("x1", -10)
    .attr("y1", function(d){
      return 0;
      //return -(1 + d.dist)/(1 + 2)*y/2;
    })
    .attr("x2", +10)
    .attr("y2", function(d){
      return 0;
      //return -(1 + d.dist)/(1 + 2)*y/2;
    })
    .attr("stroke", function (d) {
      return d.color;
    });
    
    g
    .append("polygon")
    .attr("points", "-25,0 0,25 25,0 0,-25")
    .attr("fill", function(d){
      return d.color;
    });
    
    /*
    g
    .append("polygon")
    .attr("points", function(d){
      return (-track.scale(d.size) + ",0 0,5 " + track.scale(d.size) + ",0 0,-5");
    })
    .attr("fill", function(d){
      return d.color;
    });
    */
    
    
    /*
    function (d) {
      return 0;
      //return xScale(d.start(d, i));
      //return xScale(d.mid - d.w / 2);
    })
    .attr("y", -y*0.4)
    .attr("width", function (d) {
      return d.w;
    })
    .attr("height", y*0.2)
    .attr("stroke", function (d) {
      return d.color;
    })
    .attr("fill", function (d) {
      return d.color;
    });
    */
    



  });
  
  /*
     arrow_feature.distribute (function (blocks) {
        var xScale = arrow_feature.scale();
    	blocks
    	    .select("rect")
    	    .attr("width", function (d) {
                return (xScale(d.end) - xScale(d.start));
    	    });
    });
    */
  
  // Move
  arrow_feature.move (function (arrows) {
    var xScale = arrow_feature.scale();
    var track = this;
    var y = track.height();

    arrows.select("g")
    .attr("transform", function (d) {
      return "translate(" + xScale(d.mid) + "," + (y - (0.2 + d.dist)/(1 + 2)*y) + ")";
    });
  });
  
  
  
  // Data track
  var arrow_track = tnt.board.track()
  .height(60)
  .color("white")
  .display(arrow_feature)
  .data(tnt.board.track.data.sync()
        .retriever (function () {
          return [
            {
              mid : 1,
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
  .add_track([axis_track, arrow_track]);
  
  board(div);
  
  board.start();
};
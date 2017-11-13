var custom_feature = function (div, win_size) {
  
  var board = tnt.board().from(0).to(win_size*4).min(-win_size).max(win_size*6).width(500).zoom_out(win_size*8).zoom_in(win_size);



  // Axis track
  var axis_track = tnt.board.track()
  .height(30)
  .color("white")
  .display(tnt.board.track.feature.axis()
           .orientation("bottom")
  );
  
  /*
  board.pane
  .append("text")
  .text("debug");
  */
  
  // polygon feature
  var polygon_feature = tnt.board.track.feature();

  
    
  // Create
  polygon_feature.create (function (elems) {
    
   var ymax = board.to() - board.from();
    var track = this;
    var xScale = polygon_feature.scale();
    var yScale = d3.scale.linear();
    var xresize = function(size){
      return board.width() * (size / (board.to() - board.from()));
    };
    var y = track.height();
    //var x = track.width();
    
 
    var format_xy = function(a, b){
      //console.log(xScale((a + b) / 2) + "," + (y - (2 + Math.abs(a - b))/4 * y));
      return xScale((a + b) / 2) + "," + (y - (0 + Math.abs(a - b))/ymax * y);
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
  polygon_feature.move (function (blocks) {
    var xScale = polygon_feature.scale();
    var track = this;
    var y = track.height();
    //console.log(board.to() - board.from());

    var ymax = (board.width() / (xScale(win_size) - xScale(0)) * win_size);//board.to() - board.from();
    console.log(board.width() + ' - ' +(xScale(win_size) - xScale(0)) + " - " + win_size);
    var format_xy = function(a, b){
      return xScale((a + b) / 2) + "," + (y - (0 + Math.abs(a - b))/ymax * y);
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
  

  
  
  // Data track
  var polygon_track = tnt.board.track()
  .height(300)
  .color("white")
  .display(polygon_feature)
  .data(tnt.board.track.data.sync()
        .retriever (function () {
          return [
            {
              s1 : 0,
              e1 : win_size,
              s2 : 0,
              e2 : win_size,
              color : "red",
          
            },
            {
              s1 : win_size,
              e1 : win_size*2,
              s2 : win_size,
              e2 : win_size*2,
              color : "orange",
          
            },
            {
              s1 : win_size*2,
              e1 : win_size*3,
              s2 : win_size*2,
              e2 : win_size*3,
              color : "blue",
          
            },
            {
              s1 : win_size*0,
              e1 : win_size*1,
              s2 : win_size*1,
              e2 : win_size*2,
              color : "darkred",
          
            },
            {
              s1 : win_size*1,
              e1 : win_size*2,
              s2 : win_size*2,
              e2 : win_size*3,
              color : "darkorange",
          
            },
            {
              s1 : win_size*0,
              e1 : win_size*1,
              s2 : win_size*2,
              e2 : win_size*3,
              color : "green",
          
            }
            ];
        })
  );
  
  
  board
  .add_track([polygon_track, axis_track]);
  
  board(div);
  
  board.start();
};

var custom_feature = function (div) {
  
  var board = tnt.board().from(80).to(420).min(0).max(1000).width(1000);
  
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

    var g = elems
    .append("g")
    .attr("transform", function (d) {
      return "translate(" + xScale(d.pos) + "," + y + ")";
    });
    g
    .append("line")
    .attr("x1", 0)
    .attr("y1", 0)
    .attr("x2", 0)
    .attr("y2", -y*0.9)
    .attr("stroke", function (d) {
      return d.color;
    });
    
    g
    .append("rect")
    .attr("x", function (d) {
      return -d.w/2;
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
    
    
    g
    .append("line")
    .attr("x1", 0)
    .attr("y1", 0)
    .attr("x2", -5)
    .attr("y2", -(y/2))
    .attr("stroke", function (d) {
      return d.color;
    });
    
    g
    .append("line")
    .attr("x1", 0)
    .attr("y1", 0)
    .attr("x2", 5)
    .attr("y2", -(y/2))
    .attr("stroke", function (d) {
      return d.color;
    });
  });
  
  // Move
  arrow_feature.move (function (arrows) {
    var track = this;
    var y = track.height();
    var xScale = arrow_feature.scale();
    
    arrows
    .select("g")
    .attr("transform", function (d) {
      return "translate(" + xScale(d.pos) + "," + y + ")";
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
              pos : 200,
              w : 40,
              color : "blue"
          
            },
            {
              pos : 355,
              w : 40,
              color : "orange"
            },
            {
              pos : 100,
              w : 40,
              color : "brown"
            },
            {
              pos : 400,
              w : 40,
              color : "red"
            }
            ];
        })
  );
  
  // block feature
  var block_feature = tnt.board.track.feature();
  

  
  // Create
  block_feature.create (function (elems) {
    var track = this;
        var xScale = feature.scale();
    	new_elems
    	    .append("rect")
    	    .attr("x", function (d, i) {
        		// TODO: start, end should be adjustable via the tracks API
        		return xScale(d.start(d, i));
    	    })
    	    .attr("y", 0)
    	    .attr("width", function (d, i) {
        		return (xScale(d.end(d, i)) - xScale(d.start(d, i)));
    	    })
    	    .attr("height", track.height())
    	    .attr("fill", track.color())
    	    .transition()
    	    .duration(500)
    	    .attr("fill", function (d) {
        		if (d.color === undefined) {
        		    return feature.color();
        		} else {
        		    return d.color;
        		}
    	    });
    });
  
  
      block_feature.distribute(function (elems) {
        var xScale = feature.scale();
    	elems
    	    .select("rect")
    	    .attr("width", function (d) {
        		return (xScale(d.end) - xScale(d.start));
    	    });
    });
  
  // Move
  block_feature.move (function (arrows) {
        var xScale = feature.scale();
    	blocks
    	    .select("rect")
    	    .attr("x", function (d) {
        		return xScale(d.start);
    	    })
    	    .attr("width", function (d) {
        		return (xScale(d.end) - xScale(d.start));
    	    });
  });
  
  var block_track = tnt.board.track().height(30).color("green").display(block_feature)
  .data(tnt.board.track.data.sync()
        .retriever (function () {
          return [
            {
              start : 180,
              end : 220,
              color : "blue"
          
            }
            ];
        })
  );
  
  board
  .add_track([axis_track, arrow_track, block_track]);
  
  board(div);
  
  board.start();
};
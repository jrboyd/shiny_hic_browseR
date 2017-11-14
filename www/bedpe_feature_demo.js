var bedpe_feature_demo = function (div, win_size) {
  
  var board = tnt.board().from(0).to(win_size*4).min(-win_size).max(win_size*6).width(500).zoom_out(win_size*8).zoom_in(win_size);

  // Axis track
  var axis_track = tnt.board.track()
  .height(30)
  .color("white")
  .display(tnt.board.track.feature.axis()
           .orientation("bottom")
  );
  
  var bedpe_feature = new_bedpe(board, win_size);
  
  // Data track
  var bedpe_track = tnt.board.track()
  .height(300)
  .color("white")
  .display(bedpe_feature)
  .data(tnt.board.track.data.sync()
        .retriever (function () {
          return [
            {
              s1 : 0,
              e1 : win_size*2,
              s2 : 0,
              e2 : win_size*2,
              color : "red",
          
            },
            {
              s1 : win_size*3,
              e1 : win_size*4,
              s2 : win_size,
              e2 : win_size*2,
              color : "orange",
          
            },
            {
              s1 : win_size*2,
              e1 : win_size*4,
              s2 : win_size*4,
              e2 : win_size*5,
              color : "blue",
          
            },
            {
              s1 : win_size*5,
              e1 : win_size*6,
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
  .add_track([bedpe_track, axis_track]);
  
  board(div);
  
  board.start();
};

convert                                                  \
  -delay 5                                              \
   $(for i in $(seq 0 299); do echo ${i}_frame.png; done) \
  -loop 0                                                \
   ten_ball.gif

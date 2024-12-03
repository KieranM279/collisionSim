setwd('~/Documents/collisionSim/')

library('ggplot2')

NUM_PARTICLES <- 10
GRID_SIZE <- 100
FPS <- 60
TIME <- 5
MAX_VELOCITY <- 5
RADII_RANGE <- c(1,3)
RUNTIME <- 10


temp_text <- "this is temp \n text for the\n legend of the\ngif"


plotter <- function(frame) {
  
  filename <- paste('output_data/',as.character(frame),'_frame_data.csv',sep="")
  data <- read.csv(filename)
  colnames(data)[1] <- c("Creature_ID")
  
  
  plot <- ggplot(data,aes(x=px,y=py)) +
    #xlim(0,GRID_SIZE) +
    scale_y_continuous(limits = c(0,100),breaks = seq(0,100,1)) +
    scale_x_continuous(limits = c(0,100),breaks = seq(0,100,1)) +
    geom_point(aes(size = r)) +
    theme_bw() +
    geom_hline(yintercept = 100) +
    geom_vline(xintercept = 100) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=15),
          legend.background = element_rect(fill="lightblue",
                                           linewidth=0.5, linetype="solid", 
                                           colour ="darkblue")) +
      guides(size="none")
  
  
  
  
  return(plot)
  
  
}

test <- plotter(0)

print(test)

plot_saver <- function(p,frame) {
  
  filename <- paste('frames/',as.character(frame),"_frame.png",sep = "")
  
  ggsave(plot = p,
         filename = filename,
         device = 'png',
         width = 8.7,
         height = 8.7)
}

plot_saver(p = test,frame = 0)

num_frames = FPS * TIME
num_frames = num_frames - 1

for (i in 0:num_frames) {
  print(i)
  plot <- plotter(i)
  plot_saver(p = plot,
             frame = i)
}





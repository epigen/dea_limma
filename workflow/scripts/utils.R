### Utility Functions

# extended ggsave
ggsave_new <- function(filename, results_path, plot, width=5, height=5){
    # make result directory if not exist
    if (!dir.exists(results_path)){
        dir.create(results_path, recursive = TRUE)
    }
    
    for (format in c('png')){
        ggsave(
          paste0(filename,'.',format),
          plot = plot,
          device = format,
          path = file.path(results_path),
          scale = 1,
          dpi = 300,
            width = width,
            height = height,
          limitsize = FALSE,
            units="in"
        )
    }
}


# custom theme settings
font <- "Arial"
size <- 10

custom_theme <- theme_bw(
        base_size=size,
        base_family = font
    ) %+replace% 
    
    theme(
      #text elements
      plot.title = element_text(             #title
                   family = font,            #set font family
                   size = size,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0.5,                #center align
                   vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
                   family = font,            #font family
                   size = size),               #font size
      
      plot.caption = element_text(           #caption
                   family = font,            #font family
                   size = size,                 #font size
                   hjust = 0.5),               #center align
      
      axis.title = element_text(             #axis titles
                   family = font,            #font family
                   size = size),               #font size
      
      axis.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = size),                #font size
    )
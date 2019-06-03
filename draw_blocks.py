#coding=utf-8

# import sys
import svgwrite

    ### CONSTANTS

    ## Margins
    # Bottom margin
    BOTTOM_MARGIN = 10
    # Top margin
    TOP_MARGIN = 10
    # Left margin
    LEFT_MARGIN = 30
    # Right margin
    RIGHT_MARGIN = 30

    ## Axis
    # Abcisse height from bottom margin
    ABSCISSE = 20
    # Ordonee position from the left margin
    ORDINATE = 300
    # Graduation OFFSET
    OFFSET = 10
    # Size of the graduation
    GRADSIZE = 4

    ## Mapping boxes dimensions
    # Height of the boxes
    BOX_HEIGHT = 14 
    # Space between boxes
    SPACER = 2

    ## Colors
    # BLUE
    BOX_BLUE    = svgwrite.rgb(20, 20,100, '%')
    BORDER_BLUE = svgwrite.rgb(10, 10, 66, '%')
    # RED
    BOX_RED    = svgwrite.rgb(100, 20,20, '%')
    BORDER_RED = svgwrite.rgb(66, 10, 10, '%')
    # GREEN
    BOX_GREEN    = svgwrite.rgb(20, 100 ,20, '%')
    BORDER_GREEN = svgwrite.rgb(10, 66, 10, '%')



class lane():

    def __init__(self, read_name, orientation, display_type = "block", blocks, font_size = 10, height = 14, name_space, field_sep = 5 , l_margin = 5 ):

        # read info
        self.read_name = read_name
        self.orientation = orientation
        self.blocks = blocks

        # display info
        self.display_type   = display_type
        self.font_size      = font_size
        self.height         = max(height,font_size)
        self.name_space     = name_space
        self.field_sep      = field_sep
        self.l_margin       = l_margin 
        self.graph_offset   = l_margin + name_space + field_sep
        
        # sgwrite group
        self.group = svgwrite.g(id= read_name, font_size = font_size)
        


    def get_group(self):
        return(self.group)

    def compute_blocks(self):

        pileup = {}
        for block in self.blocks:
            
            # Forward
            if(self.orientation == "0"):        
                start_b = block[0]
                end_b   = block[-1]
            # Reverse
            else:
                start_b = block[-1]
                end_b   = block[0]

            .draw_box( (offset + start_b[0], height),
                abs(end_b[0]-start_b[0]) ,
                dwg.BOX_HEIGHT,
                line_color= line_col,
                color = col,
                direction = direction,
                block_start = start_b[1],
                block_end = end_b[1]
                )

            # computing position data for this lane
            for pu in range(start_b[0],end_b[0]+1):
                pileup[pu]+=1

        return(pileup)


    def compute_dash(self):
        raise NotImplementedError("Not implemented yet")



class draw_map():

    def __init__(self, name, font_size = 10, font_color = "black", line_width = 2, tick_spacing = 25, tick_height = 4, box_spacing = 2, box_height = 15):
        
        self.dwg = svgwrite.Drawing( name+".svg" )
        self.name = name
        self.font_size  = font_size
        self.font_color = font_color
        self.line_width = line_width
        self.tick_spacing = tick_spacing
        self.tick_height = tick_height
        self.box_spacing = box_spacing
        self.box_height  = box_height
        self.origin = (0,0)

class lane():

    def __init__()


    def get_draw(self):
        return(self.dwg)


    def draw_x_axis(self, length, spacing = 25 , height = 4 ):

        dwg = self.dwg
        # number of ticks
        nb_ticks = length//spacing
        # adding group for axis
        hline = dwg.add(dwg.g(id='x_axis', stroke='black'))
        # drawing axis line
        hline.add( dwg.line( start = self.origin, end = (1 + self.origin[0] + length, self.origin[1] ), stroke = "black", stroke_width = self.line_width ))

        for y in range(nb_ticks):
            hline.add(dwg.line(start=( self.origin[0]+ (1+y)*spacing, self.origin[1] ),
                               end  =( self.origin[0]+ (1+y)*spacing, self.origin[1]-height),
                               stroke_width= self.line_width ))

            hline.add(dwg.text( str((1+y)*spacing),
                                insert=( -self.font_size + self.origin[0]+ (1+y)*spacing, self.origin[1] + self.font_size),
                                stroke = self.font_color,
                                font_size  = str(self.font_size)+'px',
                                font_family="Arial"))


    def set_size(self, width, height):
        self.dwg.size=( width, height )


    def set_bg(self,color = "white"):
        # background color
        self.dwg.add(self.dwg.rect(insert=(0, 0), size=('100%', '100%'), rx=None, ry=None, fill= color ))


    def set_origin(self,coordinates):
        self.origin = coordinates


    def draw_y_axis(self, length):

        dwg = self.dwg

        vline = dwg.add(dwg.g(id='y_axis', stroke='black'))    
        # drawing axis line
        vline.add( dwg.line( start = (self.origin[0],self.origin[1]+1), end = (self.origin[0],self.origin[1]-length), stroke = "black", stroke_width= self.line_width ))


    def draw_text(self, pos, text, text_col = "black", text_font = 'Arial' ):
        
        self.dwg.add( self.dwg.text( text,
                                insert= pos,
                                stroke = text_col,
                                font_size  = str(self.font_size)+'px',
                                font_family= text_font))


    def draw_dash(self, pos, height, line_color = "black"):

        dwg = self.dwg
        # drawing main box
        dwg.add(dwg.line(start= pos,
                         end  = (pos[0],pos[1] + height),
                         stroke_width= 1,
                         stroke = line_color))



    def draw_box(self, pos, width, height, direction = "right", line_color = "black", color = "blue", block_start = None, block_end = None, text = None):

        dwg = self.dwg
        # drawing main box
        dwg.add(dwg.rect(insert= pos,
                         size  = (width,height),
                         stroke_width= 1,
                         stroke = line_color,
                         fill= color ))
        
        # add an arrow inside the box
        if(direction):
            arrow_pos = ( pos[0] + (width // 2) , pos[1] + (height//2) )
            self.draw_arrow( arrow_pos , direction = direction, color = line_color )



        # add the relative positions of the reads
        if(block_start):
            dwg.add( dwg.text(str(block_start),
                     insert=( pos[0] + 2,  pos[1] + self.font_size ),
                     stroke = line_color,
                     font_size  = str(self.font_size)+'px',
                     font_family="Arial"))
        if(block_end):
            dwg.add( dwg.text(str(block_end),
                     insert=( pos[0] + width - len(str(block_end)) * (self.font_size - 3 ) ,  pos[1] + self.font_size ),
                     stroke = line_color,
                     font_size  = str(self.font_size)+'px',
                     font_family="Arial"))

    def draw_arrow(self, pos, offset = 5 ,direction = "right", color = "blue"):

        dwg = self.dwg

        x_offset = -offset if direction == "right" else offset

        dwg.add(dwg.line(start = pos ,
                         end   =( pos[0] + x_offset, pos[1] + x_offset ),
                         stroke_width= self.line_width,
                         stroke = color ))

        dwg.add(dwg.line(start= pos ,
                         end  =( pos[0] + x_offset, pos[1] - x_offset ),
                         stroke_width= self.line_width,
                         stroke = color ))

    def save(self):
        self.dwg.save()


# Creating draw element
# dwg = svgwrite.Drawing('test.svg', size=(500,500))


# mydwg= draw_map("test")

# mydwg.set_size(2000,500)
# mydwg.set_bg("white")
# ori = (50,250)
# mydwg.set_origin(ori)

# mydwg.draw_x_axis(2000,50,5)
# mydwg.draw_y_axis(200)


# # Ref box
# mydwg.draw_box((60,50),1800,14, line_color=BORDER_GREEN, color = BOX_GREEN, direction = "")
# # Blue box
# mydwg.draw_box((55,200),100,14, line_color=BORDER_BLUE, color = BOX_BLUE, direction = "right")
# # Red bow
# mydwg.draw_box((75,150), 1200,14, line_color=BORDER_RED, color = BOX_RED, direction = "left")
# x,y = 60, 100
# space = 10
# for i in range(20):
#     mydwg.draw_dash( (x + i* space, y ), 14 , line_color = BORDER_BLUE)


# mydwg.save()

# print("OK")
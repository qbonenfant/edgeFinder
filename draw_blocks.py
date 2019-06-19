# coding=utf-8

# import sys
import svgwrite
from collections import defaultdict as dd


class lane():

    def __init__(self, drw, read_name, ori, blocks, lane_nb, graph_offset=100, display_type="blocks", height=14, field_sep=5, l_margin=5, is_ref = False):

        # read info
        self.drawing = drw
        self.read_name = read_name
        self.orientation = "left" if ori == "1" else "right"
        self.blocks = blocks

        # display info
        self.lane_nb = lane_nb
        self.display_type = display_type
        self.font_size = drw.font_size
        self.height = max(height, drw.char_height + 2)

        self.y = drw.origin[1] + drw.SPACER + \
            drw.GRADSIZE + lane_nb * (self.height + drw.SPACER)

        self.l_margin = l_margin
        self.graph_offset = graph_offset
        self.line_col = drw.BORDER_BLUE if ori == "0" else drw.BORDER_RED
        self.color = drw.BOX_BLUE if ori == "0" else drw.BOX_RED
        if(is_ref):
            self.line_col = drw.BORDER_GREEN
            self.color = drw.BOX_GREEN
        self.pileup = dd(int)

    def compute_blocks(self):
        for block in self.blocks:

            # Forward
            if(self.orientation == "right"):
                start_b = block[0]
                end_b = block[-1]
            # Reverse
            else:
                start_b = block[-1]
                end_b = block[0]
            # Drawing the box
            self.drawing.draw_box(
                (self.drawing.origin[0] + start_b[0], self.y),
                abs(end_b[0] - start_b[0]),
                self.height,
                line_color=self.line_col,
                color=self.color,
                direction=self.orientation,
                block_start=start_b[1],
                block_end=end_b[1]
            )

            # computing position data for this lane
            for pu in range(start_b[0], end_b[0] + 1):
                self.pileup[pu] += 1

    def compute_dash(self):
        for b in self.blocks:
            # Drawing the dashs
            self.drawing.draw_dash(
                (self.drawing.origin[0] + b[0], self.y),
                self.height,
                line_color=self.line_col)
            self.pileup[b[0]] += 1

    def print_name(self):

        self.drawing.draw_text((self.drawing.LEFT_MARGIN,
                                self.y + self.height - 1),
                               self.read_name,
                               text_col=self.line_col)

    def print_lane(self):

        # drawing background rect:
        col = "lightgrey" if self.lane_nb % 2 == 0 else "white"
        self.drawing.dwg.add(self.drawing.dwg.rect(insert=(0, self.y),
                                                   size=('100%', self.height),
                                                   stroke=col,
                                                   fill=col))

        self.print_name()
        if("blocks" in self.display_type):
            self.compute_blocks()
        elif("dash" in self.display_type):
            self.compute_dash()
        else:
            raise NotImplementedError("This display method does not exists")
            print(self.display_type)

    # Comparisons
    def __lt__(self, other):
        return (self.direction <= other.direction and self.name < other.name)

    def __le__(self, other):
        return (self.direction <= other.direction and self.name <= other.name)

    def __eq__(self, other):
        return (self.direction == other.direction and self.name == other.name)

    def __ne__(self, other):
        return (self.direction != other.direction and self.name != other.name)

    def __gt__(self, other):
        return (self.direction >= other.direction and self.name > other.name)

    def __ge__(self, other):
        return (self.direction >= other.direction and self.name >= other.name)


class draw_map():

    # Margins
    # Bottom margin
    BOTTOM_MARGIN = 10
    # Top margin
    TOP_MARGIN = 10
    # Left margin
    LEFT_MARGIN = 30
    # Right margin
    RIGHT_MARGIN = 30

    # Axis
    # Abcisse height from bottom margin
    ABSCISSE = 20
    # Ordonee position from the left margin
    ORDINATE = 300
    # Graduation OFFSET
    OFFSET = 10
    # Size of the graduation
    GRADSIZE = 4

    # Mapping boxes dimensions
    # Height of the boxes
    BOX_HEIGHT = 14
    # Space between boxes
    SPACER = 4

    # Colors
    # BLUE
    BOX_BLUE = svgwrite.rgb(20, 20, 100, '%')
    BORDER_BLUE = svgwrite.rgb(10, 10, 66, '%')
    # RED
    BOX_RED = svgwrite.rgb(100, 20, 20, '%')
    BORDER_RED = svgwrite.rgb(66, 10, 10, '%')
    # GREEN
    BOX_GREEN = svgwrite.rgb(20, 100, 20, '%')
    BORDER_GREEN = svgwrite.rgb(10, 66, 10, '%')

    def __init__(self, name, font_size=10, font_color="black", line_width=2, tick_spacing=25, tick_height=4, box_spacing=2, box_height=15):

        self.dwg = svgwrite.Drawing(name + ".svg")
        self.name = name
        self.font_size = font_size
        self.font_color = font_color
        self.line_width = line_width
        self.tick_spacing = tick_spacing
        self.tick_height = tick_height
        self.box_spacing = box_spacing
        self.box_height = box_height
        self.origin = (0, 0)
        self.h_factor = float(font_size / 10)
        self.char_width, self.char_height = self.textwidth()

    def textwidth(self):
        text = "A"
        try:
            import cairo
        except Exception:
            s = len(text) * self.font_size
            return (s, s)
        surface = cairo.SVGSurface('undefined.svg', 1280, 200)
        cr = cairo.Context(surface)
        cr.select_font_face('monospace', cairo.FONT_SLANT_NORMAL,
                            cairo.FONT_WEIGHT_BOLD)
        cr.set_font_size(self.font_size)
        xb, yb, width, height, xa, ya = cr.text_extents(text)
        return (width + 4, height + 4)

    def get_draw(self):
        return(self.dwg)

    def draw_x_axis(self, length, spacing=50, height=4):

        dwg = self.dwg
        # number of ticks
        nb_ticks = length // spacing

        # adjusting to font size
        length = length * self.h_factor
        h_spacing = spacing * self.h_factor
        # adding group for axis
        hline = dwg.add(dwg.g(id='x_axis', stroke='black'))
        # drawing axis line
        hline.add(dwg.line(start=((self.origin[0] - 1), self.origin[1]),
                           end=(1 + self.origin[0] + length, self.origin[1]),
                           stroke="black",
                           stroke_width=self.line_width))

        for y in range(nb_ticks + 1):
            hline.add(dwg.line(start=(self.origin[0] + (y) * h_spacing,
                                      self.origin[1]),
                               end=(self.origin[0] + (y) *
                                    h_spacing, self.origin[1] + height),
                               stroke_width=self.line_width))
            txt = str(y * spacing)
            txt_len = self.char_width * len(txt)
            hline.add(dwg.text(txt,
                               insert=(self.origin[0] - txt_len + y * h_spacing,
                                       self.origin[1] - 2),
                               stroke=self.font_color,
                               font_size=str(self.font_size) + 'px',
                               font_family="monospace"))

    def set_bg(self, color="white"):
        # background color
        self.dwg.add(self.dwg.rect(insert=(0, 0), size=(
            '100%', '100%'), rx=None, ry=None, fill=color))

    def set_origin(self, coordinates):
        self.origin = coordinates

    def draw_y_axis(self, length):

        dwg = self.dwg

        vline = dwg.add(dwg.g(id='y_axis', stroke='black'))
        # drawing axis line
        vline.add(dwg.line(start=(self.origin[0], self.origin[1] + 1),
                           end=(self.origin[0], self.origin[1] - length),
                           stroke="black",
                           stroke_width=self.line_width))

    def draw_text(self, pos, text, text_col="black", text_font='monospace'):

        self.dwg.add(self.dwg.text(text,
                                   insert=pos,
                                   stroke=text_col,
                                   font_size=str(self.font_size) + 'px',
                                   font_family=text_font))

    def draw_dash(self, pos, height, line_color="black"):

        dwg = self.dwg
        pos = (self.origin[0] + (pos[0] - self.origin[0]) *
               self.h_factor, pos[1])
        # drawing main box
        dwg.add(dwg.line(start=pos,
                         end=(pos[0], pos[1] + height),
                         stroke_width=max(1, self.h_factor),
                         stroke=line_color))

    def draw_box(self, pos, width, height, direction="", line_color=BORDER_GREEN, color=BOX_GREEN, block_start=None, block_end=None, text=None):

        width = width * self.h_factor
        pos = (self.origin[0] + (pos[0] - self.origin[0]) *
               self.h_factor, pos[1])
        dwg = self.dwg
        # drawing main box
        dwg.add(dwg.rect(insert=pos,
                         size=(width, height),
                         stroke_width=1,
                         stroke=line_color,
                         fill=color))

        # add an arrow inside the box
        if(direction and width >= self.char_width * 2  ):
            arrow_pos = (pos[0] + (width // 2), pos[1] + height)
            self.draw_arrow(arrow_pos, direction=direction, color=line_color)

        total_txt = 0
        total_txt += len(str(block_start)) if block_start else 0
        total_txt += len(str(block_end)) if block_end else 0

        if(self.char_width * (total_txt + 3) <= width):
            # add the relative positions of the reads
            if(block_start is not None):
                dwg.add(dwg.text(str(block_start),
                                 insert=(pos[0] + 2, pos[1] + height - 1),
                                 stroke=line_color,
                                 font_size=str(self.font_size) + 'px',
                                 font_family="monospace"))
            if(block_end is not None):
                dwg.add(dwg.text(str(block_end),
                                 insert=(pos[0] + width - len(str(block_end)) *
                                         (self.char_width-2),
                                         pos[1] + height - 1),
                                 stroke=line_color,
                                 font_size=str(self.font_size) + 'px',
                                 font_family="monospace"))

    def draw_arrow(self, pos, direction="right", color="blue"):

        dwg = self.dwg

        char = ">" if direction == "right" else "<"

        dwg.add(dwg.text(char,
                         insert=pos,
                         stroke=color,
                         font_size=str(self.font_size) + 'px',
                         font_family="monospace"))

    def save(self):
        self.dwg.save()

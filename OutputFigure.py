#!/usr/bin/env python

import math
import cairo
import numpy as np
import sys
import json
import ParsedHMM
import colorsys
import collections
from scipy.stats import rankdata

class OutputFigure:
    """ Set of utilities for producing output figure

    Uses Cairo to generate a PDF figure with the name specified in the configuration.

    """

    cr = None
    padding_left = 0
    padding_top = 0
    config = None
    parsed_hmm_master = None

    def __init__(self, config):
        # Load the HMM
        self.parsed_hmm_master = ParsedHMM.ParsedHMM(config)
        # Store the config and HMM
        self.config = config
        # Set up the page
        self.padding_left = config['page']['padding_left']
        self.padding_top = config['page']['padding_top']
        surface = cairo.PDFSurface(config['output']['file_name'],
                                   self.parsed_hmm_master.length + config['page']['horizontal_padding'],
                                   config['page']['height'])
        self.cr = cairo.Context(surface)
        # Draw the master sequence
        self.draw_master_sequence()
        # Add the extra hits
        self.add_hits()
        # Save the file
        self.save_file()

    def plot_clustal(self, scale, max_bitscore, position, offset_horizontal, draw_full_rectangle, hmm):
        # Draw axes
        self.cr.set_source_rgb(0, 0, 0)
        self.cr.set_line_width(1)
        self.cr.move_to(offset_horizontal - 0.5, position - 1)
        self.cr.line_to(offset_horizontal - 0.5, position + scale * max_bitscore)
        self.cr.move_to(offset_horizontal - 1, position + + scale * max_bitscore)
        self.cr.line_to(offset_horizontal + len(hmm.probs) + 0.5, position + scale * max_bitscore)
        if draw_full_rectangle:
            self.cr.move_to(offset_horizontal - 1, position - 1)
            self.cr.line_to(offset_horizontal + len(hmm.probs) + 0.5, position - 1)
            self.cr.move_to(offset_horizontal + len(hmm.probs) + 0.5, position - 1)
            self.cr.line_to(offset_horizontal + len(hmm.probs) + 0.5, position + scale * max_bitscore + 0.5)
        self.cr.set_font_size(4 * np.log(scale * max_bitscore / 5))
        self.cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        for i in range(6):
            (x, y, width, height, dx, dy) = self.cr.text_extents(str(i * max_bitscore / 5))
            self.cr.move_to(offset_horizontal - 10 - width,
                            position + height / 2 + scale * max_bitscore - (max_bitscore / 5 * i * scale))
            self.cr.show_text(str(i * max_bitscore / 5))
            self.cr.move_to(offset_horizontal - 1, position + scale * max_bitscore - (max_bitscore / 5 * i * scale))
            self.cr.line_to(offset_horizontal - 1 - 6, position + scale * max_bitscore - (max_bitscore / 5 * i * scale))

        self.cr.stroke()
        self.cr.close_path()
        self.cr.set_font_size(4 * np.log(scale * max_bitscore / 5))
        (x, y, width, height, dx, dy) = self.cr.text_extents("Information (bits)")
        self.cr.move_to(offset_horizontal - 50 - height / 2, position + (scale * max_bitscore - width) / 2)
        self.cr.rotate(math.pi / 2)
        self.cr.show_text("Information (bits)")
        self.cr.rotate(-math.pi / 2)

        # Plot the clustal plot on top
        for i in range(len(hmm.probs)):
            height_offset = 0
            for j in range(len(self.config['colours'])):
                # Calculate colour
                if np.sum(hmm.clustal_colours[i]) > 1.5:
                    col = list(colorsys.rgb_to_hsv(self.config['colours'][j]['rgb'][0], self.config['colours'][j]['rgb'][1],
                                           self.config['colours'][j]['rgb'][2]))
                    col[1] = (np.divide(np.sum(hmm.clustal_colours[i]), (np.max(hmm.height_array))))
                else:
                    col = [0, 0, 0.8]

                col = colorsys.hsv_to_rgb(col[0], col[1], col[2])
                self.cr.set_source_rgb(col[0], col[1], col[2])

                height = scale * hmm.clustal_colours[i][j]
                self.cr.rectangle(offset_horizontal + i, position - height + scale * max_bitscore + height_offset, 1,
                                  height)
                height_offset -= height
                self.cr.fill()

    def draw_master_sequence(self):
        self.plot_clustal(30/7 * 5, 7, self.padding_top, self.padding_left, self.config['output']['split'],
                          self.parsed_hmm_master)

        if self.config['output']['subplot_type'] == "secondary":
            self.draw_ss(self.parsed_hmm_master, self.padding_top + 225, self.padding_left, -1, self.parsed_hmm_master.length)
        elif self.config['output']['subplot_type'] == "psiplot":
            self.bar_ss(self.parsed_hmm_master, self.padding_top + 225, self.padding_left, -1,
                         self.parsed_hmm_master.length)

        # Draw a box to represent the main sequence
        draw_masters_at = [self.padding_top]
        if self.config['output']['split']:
            draw_masters_at.append(self.padding_top - 250)

        for vertical_pos_index in range(len(draw_masters_at)):
            vertical_pos = draw_masters_at[vertical_pos_index]
            above_or_below = np.remainder(vertical_pos_index, 2) * 1

            # Draw some boxes for KKT4
            self.cr.set_source_rgb(0.2, 0.2, 0.2)
            self.cr.set_line_width(1)
            self.cr.rectangle(self.padding_left, vertical_pos + 180, self.parsed_hmm_master.length, 40)
            self.cr.stroke()
            self.cr.close_path()
            self.cr.set_source_rgba(0, 1, 0, 0.1)
            self.cr.rectangle(self.padding_left + above_or_below * self.config['output']['split_at'],
                              vertical_pos + 180,
                              (-(above_or_below - 1)) * self.config['output']['split_at'] + above_or_below * (
                                      self.parsed_hmm_master.length - self.config['output']['split_at']), 40)
            self.cr.fill()
            self.cr.close_path()
            self.cr.set_font_size(20)
            self.cr.set_source_rgb(0, 0, 0)
            (x, y, width, height, dx, dy) = self.cr.text_extents(self.parsed_hmm_master.name)
            self.cr.move_to(self.padding_left + self.parsed_hmm_master.length + 15, vertical_pos + 200 + (height) / 2)
            self.cr.show_text(self.parsed_hmm_master.name)
            # Add on some residue numbering above the box
            self.cr.set_font_size(12)
            for i in range(0, self.parsed_hmm_master.length - 50, 100):
                (x, y, width, height, dx, dy) = self.cr.text_extents(str(i))
                self.cr.move_to(self.padding_left + i - width / 2,
                                vertical_pos + 165 + height / 2 + above_or_below * (180 - 151 + 40))
                self.cr.show_text(str(i))
                self.cr.move_to(self.padding_left + i, vertical_pos + 151 + above_or_below * (180 - 151 + 40))
                self.cr.line_to(self.padding_left + i, vertical_pos + 151 + 6 + above_or_below * (180 - 151 + 40))
                self.cr.move_to(self.padding_left + i, vertical_pos + 165 + height + above_or_below * (180 - 151 + 40))
                self.cr.line_to(self.padding_left + i,
                                vertical_pos + 165 + height + 6 + above_or_below * (180 - 151 + 40))

            (x, y, width, height, dx, dy) = self.cr.text_extents(str(self.parsed_hmm_master.length))
            self.cr.move_to(self.padding_left + self.parsed_hmm_master.length - width / 2,
                            vertical_pos + 165 + height / 2 + above_or_below * (180 - 151 + 40))
            self.cr.show_text(str(self.parsed_hmm_master.length))
            self.cr.move_to(self.padding_left + self.parsed_hmm_master.length,
                            vertical_pos + 151 + above_or_below * (180 - 151 + 40))
            self.cr.line_to(self.padding_left + self.parsed_hmm_master.length,
                            vertical_pos + 151 + 6 + above_or_below * (180 - 151 + 40))
            self.cr.move_to(self.padding_left + self.parsed_hmm_master.length,
                            vertical_pos + 165 + height + above_or_below * (180 - 151 + 40))
            self.cr.line_to(self.padding_left + self.parsed_hmm_master.length,
                            vertical_pos + 165 + height + 6 + above_or_below * (180 - 151 + 40))
            self.cr.stroke()

            # Add domains from pfam folder
            with open(self.config['searches']['pfam'], 'r') as file:
                lines = file.readlines()
                # Look for start and end position of each hit
                for hit_group in range(len(self.config['domains'])):
                    # Store start and end position
                    start = 100000
                    end = 0
                    # Move on to
                    if isinstance(self.config['domains'][hit_group]['pfam_hit_number'], list):
                        for hit in range(len(self.config['domains'][hit_group]['pfam_hit_number'])):
                            for line in lines:
                                # Look for hits of our numbers
                                if line.strip().startswith(
                                        str(self.config['domains'][hit_group]['pfam_hit_number'][hit]) + " "):
                                    columns = line.split()
                                    start = np.min([start, int(columns[len(columns) - 3].split('-')[0])])
                                    end = np.max([end, int(columns[len(columns) - 3].split('-')[1])])
                    else:
                        # Look for hits of our numbers
                        for line in lines:
                            if line.strip().startswith(str(self.config['domains'][hit_group]['pfam_hit_number']) + " "):
                                columns = line.split()
                                start = int(columns[len(columns) - 3].split('-')[0])
                                end = int(columns[len(columns) - 3].split('-')[1])

                    # Now draw that hit group on
                    if start > self.config['output']['split_at'] and vertical_pos_index is 0 or start < \
                            self.config['output'][
                                'split_at'] and vertical_pos_index is 1:
                        self.cr.set_source_rgba(0.5, 0.5, 0.5, 0.5)
                    else:
                        self.cr.set_source_rgba(self.config['domains'][hit_group]['colour'][0],
                                                self.config['domains'][hit_group]['colour'][1],
                                                self.config['domains'][hit_group]['colour'][2],
                                                self.config['domains'][hit_group]['colour'][3])
                    self.cr.set_line_width(1)
                    self.cr.rectangle(self.padding_left + start, vertical_pos + 180.5, end - start, 39)
                    self.cr.fill()
                    self.cr.stroke()
                    # Label it
                    self.cr.set_font_size(20)
                    self.cr.set_source_rgb(0, 0, 0)
                    name = self.config['domains'][hit_group]['name']
                    (x, y, width, height, dx, dy) = self.cr.text_extents(name)
                    self.cr.move_to(self.padding_left + start + ((end - start) - width) / 2,
                                    vertical_pos + 200 + height / 2)
                    self.cr.show_text(name)

    def add_hits(self):
        # Calculate a spacing factor
        if self.config['output']['subplot_type'] == "logo":
            spacing = 150
            top_offset = -1.4
        else:
            spacing = 120
            top_offset = -2.2

        e_value_array = []
        with open(self.config['searches']['hits'], 'r') as file:
            lines = file.readlines()
            # Look for start and end position of each hit
            for hit_group in range(2, self.config['output']['max_hits'] + 1):
                for line in lines:
                    if line.strip().startswith(str(hit_group) + " "):
                        columns = line.split()
                        e = float(columns[len(columns) - 8])
                        e_value_array.append(e)

        # Start adding hits on
        with open(self.config['searches']['hits'], 'r') as file:
            idx = [0.5, top_offset]  # first is below, second is above
            lines = file.readlines()

            # Store e value in correct scope to be able to compare
            last_e = 0

            # Look for start and end position of each hit
            for counter in range(2, self.config['output']['max_hits'] + 1):
                print("--- Starting new OG ---")
                print(e_value_array)
                print(rankdata(np.array(e_value_array), 'ordinal'))
                hit_group_in_file = int(np.where(rankdata(np.array(e_value_array), 'ordinal') == counter - 1)[0] + 2)
                print(hit_group_in_file)


                for line in lines:
                    if line.strip().startswith(str(hit_group_in_file) + " "):
                        columns = line.split()
                        start = int(columns[len(columns) - 3].split('-')[0])
                        end = int(columns[len(columns) - 3].split('-')[1])
                        e = float(columns[len(columns) - 8])

                        hit_start = int(columns[len(columns) - 2].split('-')[0])
                        hit_end = int(columns[len(columns) - 2].split('-')[1])
                        hit_length = int(columns[len(columns) - 1].replace("(", "").replace(")", ""))
                        name = columns[1]
                        self.cr.set_line_width(1)
                        current_idx = counter - 0.5
                        if self.config['output']['split']:
                            if start < self.config['output']['split_at']:
                                idx[0] += 1
                                current_idx = idx[0]
                            else:
                                idx[1] -= 1
                                current_idx = idx[1]

                        # Should we draw a line?
                        if not self.config['output']['split']:
                            message = ''
                            if e > 0.00001 and last_e < 0.00001:
                                message = "E = 0.00001"
                            if e >= 0.001 and last_e < 0.001:
                                message = "E = 0.001"
                            if e >= 0.05 and last_e < 0.05:
                                message = "E = 0.05"
                            if e >= 1 and last_e < 1:
                                message = "E = 1"
                            if e >= 10 and last_e < 10:
                                message = "E = 10"

                            if message is not "":
                                print(message)
                                self.cr.set_source_rgb(0, 0, 0)
                                self.cr.set_dash([30, 10, 10, 10])
                                self.cr.move_to(0, self.padding_top + 123 + spacing * current_idx)
                                self.cr.line_to(self.config['page']['padding_left'] + self.parsed_hmm_master.length +
                                                self.config['page']['horizontal_padding'], self.padding_top + 123 + spacing * current_idx)
                                self.cr.stroke()
                                self.cr.set_font_size(22)
                                (x, y, width, height, dx, dy) = self.cr.text_extents(message)
                                self.cr.move_to(self.parsed_hmm_master.length + self.config['page']['horizontal_padding'] - width - 15, self.padding_top + 123 + spacing * current_idx - height - 2)
                                self.cr.show_text(message)
                            last_e = e

                        self.cr.set_dash([])
                        self.cr.rectangle(self.padding_left + start, self.padding_top + 150 + spacing * current_idx,
                                          hit_end - hit_start, 40)
                        # If there's a cutoff, go for cutoff based e-value display
                        print('cutoff' in self.config['output'])
                        if 'cutoff' in self.config['output']:
                            if e < self.config['output']['cutoff']:
                                self.cr.set_source_rgb(0.6, 1, 0.6)
                            else:
                                self.cr.set_source_rgb(0.9, 0.9, 0.9)
                        else:
                            print("E: " + str(e))
                            max_e = np.max(e_value_array)
                            min_e = np.min(e_value_array)
                            e_range = max_e - min_e
                            col = [0, 0, 0]

                            # For mixing based on e-value, choose this
                            #col[0] = (120/360)-(120/360)*(e-min_e)/np.log(e_range)
                            col[0] = (120/360) * (1 - (counter - 2) / (len(e_value_array) - 1))

                            print("H: " + str(col[0]))
                            col[1] = 0.5
                            col[2] = 1
                            col = colorsys.hsv_to_rgb(col[0], col[1], col[2])
                            self.cr.set_source_rgb(col[0], col[1], col[2])

                        self.cr.fill()
                        self.cr.rectangle(self.padding_left + start, self.padding_top + 150 + spacing * current_idx,
                                          hit_end - hit_start, 40)
                        # Main hit region drawn

                        ## This will be useful in the future
                        og_ends_at = self.padding_left + start - hit_start + hit_length


                        # Now draw the bits that didn't hit in grey
                        self.cr.set_source_rgb(0.6, 0.6, 0.6)
                        # Check if we go near the end of page
                        if og_ends_at + 200 > self.config['page']['padding_left'] + self.parsed_hmm_master.length + self.config['page']['horizontal_padding']:
                            truncated = 400 + (og_ends_at - (self.config['page']['padding_left'] + self.parsed_hmm_master.length + self.config['page']['horizontal_padding']))
                            og_ends_at = self.config['page']['padding_left'] + self.parsed_hmm_master.length + self.config['page']['horizontal_padding'] - 400
                            self.cr.rectangle(self.padding_left + start - hit_start,
                                              self.padding_top + 150 + spacing * current_idx,
                                              og_ends_at - (self.padding_left + start - hit_start), 40)
                            self.cr.stroke()
                            # Draw an arrow and indicate how much was truncated in small writing
                            arrow_length = 60
                            arrow_angle = 0
                            arrowhead_angle = math.pi / 6
                            arrowhead_length = 20
                            self.cr.move_to(og_ends_at, self.padding_top + 150 + spacing * current_idx + 20)  # move to center of canvas

                            self.cr.rel_line_to(arrow_length * math.cos(arrow_angle), arrow_length * math.sin(arrow_angle))
                            self.cr.rel_move_to(-arrowhead_length * math.cos(arrow_angle - arrowhead_angle),
                                            -arrowhead_length * math.sin(arrow_angle - arrowhead_angle))
                            self.cr.rel_line_to(arrowhead_length * math.cos(arrow_angle - arrowhead_angle),
                                            arrowhead_length * math.sin(arrow_angle - arrowhead_angle))
                            self.cr.rel_line_to(-arrowhead_length * math.cos(arrow_angle + arrowhead_angle),
                                            -arrowhead_length * math.sin(arrow_angle + arrowhead_angle))

                            self.cr.set_source_rgb(0, 0, 0)
                            self.cr.set_line_width(2)
                            self.cr.stroke()

                            continuation_text = "+" + str(truncated) + "aa"
                            self.cr.set_font_size(18)
                            (x, y, width, height, dx, dy) = self.cr.text_extents(continuation_text)
                            self.cr.move_to(og_ends_at + 5,
                                            self.padding_top + 150 + spacing * current_idx - height + 12)


                            self.cr.show_text(continuation_text)

                            # Update new end position
                            og_ends_at = og_ends_at + 70
                        else:
                            # Not near end of page
                            self.cr.rectangle(self.padding_left + start - hit_start, self.padding_top + 150 + spacing * current_idx,
                                              hit_length,
                                              40)
                        self.cr.stroke()
                        # Draw a name on
                        self.cr.set_source_rgb(0, 0, 0)
                        self.cr.set_font_size(20)
                        (x, y, width, height, dx, dy) = self.cr.text_extents(name)
                        if width < hit_end-hit_start:
                            # Place OG name within hit region
                            self.cr.move_to(self.padding_left + start + ((hit_end - hit_start) - width) / 2,
                                            self.padding_top + 150 + 20 + spacing * current_idx + height / 2)
                        elif width < hit_start:
                            # Hit region box starts at self.padding_left + start
                            # Overall hit bix starts at self.padding_left + start - hit_start

                            # Can place in leftmost white box?
                            self.cr.move_to(self.padding_left + start - hit_start - width / 2 + hit_start/2,
                                            self.padding_top + 150 + 20 + spacing * current_idx + height / 2)
                        elif (width < hit_length - hit_end):
                            # Can place in rightmost white box?
                            self.cr.move_to(((og_ends_at + (self.padding_left + start + hit_end - hit_start )) - width) / 2,
                                            self.padding_top + 150 + 20 + spacing * current_idx + height / 2)
                        elif width < (self.padding_left + start):
                            print("Place left")
                            # Place on left of box
                            self.cr.move_to(self.padding_left + start - hit_start - width - 5, self.padding_top + 150 + 20 + spacing * current_idx + height / 2)
                        else:
                            print("Place right")
                            # Place on right of box
                            self.cr.move_to(self.padding_left + start - hit_start+hit_length, self.padding_top + 150 + 20 + spacing * current_idx + height / 2)

                        self.cr.show_text(name)
                        self.cr.stroke()
                        # Load the HMM
                        hmm = ParsedHMM.ParsedHMM(
                            {"master": {"name": name, "hmm_file": "hmms/" + name + ".fa.hmm.ss.hmm"},
                             "colours": self.config['colours']})
                        print(name)
                        print(hmm.length)
                        # Plot either logo or 2ndary structure
                        if self.config['output']['subplot_type'] == "logo":
                            max_height = np.ceil(np.average(hmm.height_array) / 10) * 10
                            scale = (30 / 5) * (10 / max_height)
                            # Add the clustal plot
                            self.plot_clustal(scale, max_height, self.padding_top + 90 + spacing * current_idx,
                                              self.padding_left + start - hit_start, False, hmm)
                        elif self.config['output']['subplot_type'] == "secondary":
                            self.draw_ss(hmm, self.padding_top + 195 + spacing * current_idx, self.padding_left + start - hit_start, hit_start, hit_end, og_ends_at)
                        elif self.config['output']['subplot_type'] == "psiplot":
                            self.bar_ss(hmm, self.padding_top + 195 + spacing * current_idx, self.padding_left + start - hit_start, hit_start, hit_end, og_ends_at)



    def bar_ss(self, hmm, pos_y, pos_x, hit_start, hit_end, right_cutoff=100000):
        # (self, scale, max_bitscore, position, offset_horizontal, draw_full_rectangle, hmm):
        # H is green, C is yellow, E is grey
        max_bitscore = 9
        offset_horizontal = pos_x
        position = pos_y
        print(list(set(hmm.ss)))
        scale = 3
        # Draw axes
        self.cr.set_source_rgb(0, 0, 0)
        self.cr.set_line_width(1)
        self.cr.move_to(offset_horizontal - 0.5, position - 1)
        self.cr.line_to(offset_horizontal - 0.5, position + scale * max_bitscore)
        self.cr.move_to(offset_horizontal - 1, position + + scale * max_bitscore)
        self.cr.line_to(offset_horizontal + len(hmm.probs) + 0.5, position + scale * max_bitscore)
        if True:
            self.cr.move_to(offset_horizontal - 1, position - 1)
            self.cr.line_to(offset_horizontal + len(hmm.probs) + 0.5, position - 1)
            self.cr.move_to(offset_horizontal + len(hmm.probs) + 0.5, position - 1)
            self.cr.line_to(offset_horizontal + len(hmm.probs) + 0.5, position + scale * max_bitscore + 0.5)
        self.cr.set_font_size(4 * np.log(scale * max_bitscore / 2))
        self.cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        for i in range(3):
            (x, y, width, height, dx, dy) = self.cr.text_extents(str(i * max_bitscore/2))
            self.cr.move_to(offset_horizontal - 10 - width,
                            position + height / 2 + scale * max_bitscore - (max_bitscore / 2 * i * scale))
            self.cr.show_text(str(i * max_bitscore/2))
            self.cr.move_to(offset_horizontal - 1, position + scale * max_bitscore - (max_bitscore / 2 * i * scale))
            self.cr.line_to(offset_horizontal - 1 - 6, position + scale * max_bitscore - (max_bitscore / 2 * i * scale))

        self.cr.stroke()
        self.cr.close_path()
        self.cr.set_font_size(4 * np.log(scale * max_bitscore / 2))
        (x, y, width, height, dx, dy) = self.cr.text_extents("SS Conf.")
        self.cr.move_to(offset_horizontal - 50 - height / 2, position + (scale * max_bitscore - width) / 2)
        self.cr.rotate(math.pi / 2)
        self.cr.show_text("SS Conf.")
        self.cr.rotate(-math.pi / 2)

        print(hmm.ss_probs)

        # Plot the clustal plot on top
        for i in range(len(hmm.probs)):
            height_offset = 0

            #helix SS is in red. The sheet SS is in green. The coil SS is in gray
            col = [0, 1, 0]
            if hmm.ss[i] == "H":
                col = [1, 0, 0]
            elif hmm.ss[i] == "C":
                col = [0.5, 0.5, 0.5]

            # Height is the ss_prob
            height = int(hmm.ss_probs[i]) * scale

            self.cr.set_source_rgb(col[0], col[1], col[2])

            self.cr.rectangle(offset_horizontal + i, position - height + scale * max_bitscore + height_offset, 1,
                              height)
            height_offset -= height
            self.cr.fill()


    def draw_ss(self, hmm, pos_y, pos_x, hit_start, hit_end, right_cutoff=100000):
        # First of all, calculate what moving average we'll need to take
        font_size = 12
        self.cr.set_font_size(font_size)
        (x, y, width, height, dx, dy) = self.cr.text_extents("W")
        pos_y = pos_y + height
        bin_size = math.floor(width)

        # Go through in each bin size
        idx = 0
        while idx < hmm.length/bin_size:
            bin_mean_sequence_pos = idx*bin_size - bin_size/2

            ss_in_bin = hmm.ss[idx*bin_size:(idx+1)*bin_size]
            # print(ss_in_bin)
            mean_structure = list(collections.Counter(ss_in_bin).most_common(1)[0])[0]
            confidence_in_bin = hmm.ss_probs[idx*bin_size:(idx+1)*bin_size]
            # print(confidence_in_bin)
            sum_confidence = 0
            for c in list(confidence_in_bin):
                sum_confidence += int(c)
            bin_confidence = int(np.round(sum_confidence/bin_size))

            if pos_x + bin_mean_sequence_pos > right_cutoff:
                return

            # Now use mean_structure, bin_confidence to plot
            if bin_mean_sequence_pos < hit_start or bin_mean_sequence_pos > hit_end:
                self.cr.set_source_rgb(0.2, 0.2, 0.2)
            else:
                self.cr.set_source_rgb(0, 0, 1)

            self.cr.move_to(pos_x + bin_mean_sequence_pos,
                            pos_y)
            self.cr.show_text(mean_structure)

            # If not in the hit range, show confidence in grey
            if bin_mean_sequence_pos < hit_start or bin_mean_sequence_pos > hit_end:
                col = [0, 0, 1]
                col[2] = 1-np.divide(bin_confidence, 9)
            else:
                # If inside hit range, confidence in green
                col = [0, 0, 1]
                col = list(colorsys.rgb_to_hsv(col[0], col[1], col[2]))
                # Confidence goes from 0 to 9
                col[1] = np.divide(bin_confidence, 9)

            col = colorsys.hsv_to_rgb(col[0], col[1], col[2])
            self.cr.set_source_rgba(col[0], col[1], col[2])
            self.cr.rectangle(pos_x + bin_mean_sequence_pos, pos_y+3, bin_size, height+2)
            self.cr.fill()
            # Draw white number in
            self.cr.set_source_rgba(1, 1, 1, 1)
            self.cr.move_to(pos_x + bin_mean_sequence_pos + 2, pos_y+height+4)
            self.cr.show_text(str(bin_confidence))

            idx += 1

    def save_file(self):
        self.cr.save()
        self.cr.show_page()

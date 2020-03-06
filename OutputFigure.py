#!/usr/bin/env python

import math
import cairo
import numpy as np
import sys
import json
import ParsedHMM


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
                if np.sum(hmm.clustal_colours[i]) > np.percentile(
                        hmm.height_array, 5):
                    self.cr.set_source_rgb(self.config['colours'][j]['rgb'][0], self.config['colours'][j]['rgb'][1],
                                           self.config['colours'][j]['rgb'][2])
                else:
                    self.cr.set_source_rgb(0.3, 0.3, 0.3)
                height = scale * hmm.clustal_colours[i][j]
                self.cr.rectangle(offset_horizontal + i, position - height + scale * max_bitscore + height_offset, 1,
                                  height)
                height_offset -= height
                self.cr.fill()

    def draw_master_sequence(self):
        self.plot_clustal(30, 5, self.padding_top, self.padding_left, self.config['output']['split'],
                          self.parsed_hmm_master)

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
        # Start adding hits on
        with open(self.config['searches']['hits'], 'r') as file:
            idx = [0.5, -1.4]  # first is below, second is above
            lines = file.readlines()
            # Look for start and end position of each hit
            for hit_group in range(2, self.config['output']['max_hits'] + 1):
                for line in lines:
                    if line.strip().startswith(str(hit_group) + " "):
                        columns = line.split()
                        start = int(columns[len(columns) - 3].split('-')[0])
                        end = int(columns[len(columns) - 3].split('-')[1])
                        e = float(columns[len(columns) - 8])
                        hit_start = int(columns[len(columns) - 2].split('-')[0])
                        hit_end = int(columns[len(columns) - 2].split('-')[1])
                        hit_length = int(columns[len(columns) - 1].replace("(", "").replace(")", ""))
                        name = columns[1]
                        self.cr.set_line_width(1)
                        current_idx = hit_group - 0.5
                        if self.config['output']['split']:
                            if start < self.config['output']['split_at']:
                                idx[0] += 1
                                current_idx = idx[0]
                            else:
                                idx[1] -= 1
                                current_idx = idx[1]

                        self.cr.rectangle(self.padding_left + start, self.padding_top + 150 + 150 * current_idx,
                                          np.min([end - start, hit_length - hit_start]), 40)
                        if e < self.config['output']['cutoff']:
                            self.cr.set_source_rgb(0.6, 1, 0.6)
                        else:
                            self.cr.set_source_rgb(0.9, 0.9, 0.9)
                        self.cr.fill()
                        self.cr.rectangle(self.padding_left + start, self.padding_top + 150 + 150 * current_idx,
                                          np.min([end - start, hit_length - hit_start]), 40)
                        # Main hit region drawn
                        # Now draw the bits that didn't hit in grey
                        self.cr.set_source_rgb(0.6, 0.6, 0.6)
                        self.cr.rectangle(self.padding_left + start - hit_start, self.padding_top + 150 + 150 * current_idx,
                                          hit_length,
                                          40)
                        self.cr.stroke()
                        # Draw a name on
                        self.cr.set_source_rgb(0, 0, 0)
                        self.cr.set_font_size(20)
                        (x, y, width, height, dx, dy) = self.cr.text_extents(name)
                        self.cr.move_to(self.padding_left + start + ((end - start) - width) / 2,
                                        self.padding_top + 150 + 20 + 150 * current_idx + height / 2)
                        self.cr.show_text(name)
                        self.cr.stroke()
                        # Load the HMM
                        hmm = ParsedHMM.ParsedHMM(
                            {"master": {"name": name, "hmm_file": "hmms/" + name + ".fa.hmm.ss.hmm"},
                             "colours": self.config['colours']})
                        print(name)
                        print(hmm.length)
                        max_height = np.ceil(np.average(hmm.height_array) / 10) * 10
                        scale = (30 / 5) * (10 / max_height)
                        # Add the clustal plot
                        self.plot_clustal(scale, max_height, self.padding_top + 90 + 150 * current_idx,
                                          self.padding_left + start - hit_start, False, hmm)

    def save_file(self):
        self.cr.save()
        self.cr.show_page()

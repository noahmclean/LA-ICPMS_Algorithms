#!/usr/bin/env python2.7
# Copyright (c) 2014 Dr. Philip Wenig
# Copyright (c) 2015 John H. Hartman
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License version
# 2.1, as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License version 2.1 along with this program.
# If not, see <http://www.gnu.org/licenses/>.
#
import struct
from pprint import *
import sys
import os
import optparse
import glob
import datetime
import math

HDR_INDEX_OFFSET = 33
HDR_INDEX_LEN = 39
HDR_TIMESTAMP = 40

SCAN_NUMBER = 9
SCAN_DELTA = 7
SCAN_ACF = 12
SCAN_PREV_TIME = 18
SCAN_TIME = 19

scales = {}

for i in xrange(0,16):
    key = 0x1010 + i
    scales[key] = math.pow(2, i) * 1.0
    #print "scales[0x%x] = %f" % (key, scales[key])

for scale in scales.keys():
    key = scale & 0xFF0F
    #print "key = 0x%x" % key
    scales[key] = scales[scale]
    #key = scale | 0x400
    #scales[key] = scales[scale]

def Scale(value, scale, acf):
    #print "Scale(%d, 0x%x, %f)" % (value, scale, acf)
    if scale & 0xF00:
        # Don't know what this digit does, flag the result
        suffix = '*'
        scale = scale & 0xF0FF
    else:
        suffix = ''
    if scale not in scales:
        raise Exception("Unknown scaling %s (%d)" % (hex(scale), value))
    factor = scales[scale]
    if scale & 0XF0 == 0:
        factor *= acf
    result = str(value * factor) + suffix
    #print "result = %s" % result
    return result

magic = 0x8000

def Magic(value):
    global magic
    if value == magic:
        return True
    elif value == magic + 1:
        magic += 1
        return True
    else:
        return False

def ReadHeader(dat):

    dat.seek(0x10)
    fields = 85
    data = dat.read(fields * 4)
    hdr = struct.unpack('<%dI' % fields, data)
    return hdr

def main(args):

    usage = "Usage %prog [options] [input files]"

    parser = optparse.OptionParser(version="%prog 1.4", usage=usage)
    parser.add_option("-c", "--comments",
                      action="store_true", dest="comments",
                      default=False,
                      help="add diagnostic comments to the output file")
    (options, args) = parser.parse_args(args[1:])

    parser.print_version()
    inputs = []
    if len(args) == 0:
        args = ["untSMPABC001.dat"]

    # Expand any directories into the dat files they contain.

    dirs = []
    for arg in args:
        if os.path.isdir(arg):
            dirs.append(arg)
            inputs += glob.glob(os.path.join(arg, '*.dat'))
        else:
            inputs.append(arg)

    # Sort the dat files by their creation time. 

    timestamps = {}
    hdrs = {}
    for datfile in inputs:
        dat = open(datfile, 'rb')
        hdr = ReadHeader(dat)
        dat.close()
        timestamps[hdr[HDR_TIMESTAMP]] = datfile
        hdrs[datfile] = hdr

    inputs = [timestamps[x] for x in sorted(timestamps.keys())]

    # Determine the output directory.

    outputdir = os.path.expanduser('~/Desktop')
    if len(dirs) == 1:
        outputdir = dirs[0]
    elif len(dirs) == 0:
        # No directory specified, try to infer it from the input files.
        dirs = {}
        for datfile in inputs:
            d = os.path.split(datfile)[0]
            dirs[d] = 1
        if len(dirs) == 1:
            outputdir = dirs.keys()[0]
    base = os.path.splitext(os.path.split(inputs[0])[1])[0]
    if len(inputs) > 1:
        base += 'combined'
        i = 0
        while True:
            outputfile = os.path.join(outputdir, base + '%02d' % i + '.csv')
            if not os.path.exists(outputfile):
                break
            i += 1
    else:
        outputfile = os.path.join(outputdir, base + '.csv')
    output = open(outputfile, 'w')

    print "Writing to", outputfile
    print "Using default elements"
    elements = [202, 204, 206, 207, 208, 232, 238]
    first = True
    for datfile in inputs:

        dat = open(datfile, 'rb')
        #inf = open(base + '.inf', 'rb')

        hdr = hdrs[datfile]

        timestamp = hdr[HDR_TIMESTAMP]
        msg = ' '.join([datfile, str(timestamp), str(datetime.datetime.fromtimestamp(timestamp))])
        print msg
        if options.comments:
            print >> output, "# ", msg

        indexOffset = hdr[HDR_INDEX_OFFSET] + 4
        indexLength = hdr[HDR_INDEX_LEN]
        dat.seek(indexOffset)
        offsets = struct.unpack('<' + indexLength * 'I', dat.read(indexLength * 4))
        scans = len(offsets)
        prev = 0
        for offset in offsets:
            #print hex(offset), str(offset - prev)
            if prev != 0:
                indexSize = offset - prev
            prev = offset

        count = (indexSize - (22 * 4)) / 2

        # Read elements
        """
        inf.seek(0x32D)
        offset = struct.unpack('<H', inf.read(2))[0]
        print "offset = %d" % offset
        inf.seek(offset)
        for scan in xrange(scans):
            value = struct.unpack('<2d', inf.read(16))[0]
            print "value = %f" % value
            if (value >= 1) and (value <= 300):
                element = int(round(value, 0))
                if element not in elements:
                    print "element %d" % element
                    elements.append(element)
            else:
                break

        sys.exit(0)
        """
        for offset in offsets:
            #print "offset = %d (0x%x)" % (offset, offset)
            dat.seek(offset)
            data = dat.read(indexSize)
            vals = struct.unpack('<22I%uH' % count, data)
            #print "# vals = %d" % len(vals)
            scanNumber = vals[SCAN_NUMBER]
            #print "scan = %d" % scanNumber
            acf = (vals[SCAN_ACF] * 1.0) / 64.0
            #print "acf = %s" % acf
            #acquisitions = vals[8]
            tmp = [str(vals[SCAN_TIME]/ 1000.0)]
            base = 72
            key = vals[base:base+2]
            #print "key = " + str(key)
            t = vals[SCAN_TIME]/ 1000.0 + timestamp
            #print "time = %f" % t
            #result = [str(scanNumber), str(t), str(vals[SCAN_DELTA] / 1000.0), str(vals[SCAN_PREV_TIME] / 1000.0), str(acf)]
            #headers = ["Scan", "Time", "Delta?", "Prev Time?", "ACF"]
            result = [str(scanNumber), '%f' % t, '%f' % acf]
            headers = ["Scan", "Time", "ACF"]
            pulses = []
            analogs = []
            mysteries = ['']*5
            averages = [str(t)]
            index = base
            total = 0.0
            n = 0
            acquisition = 0
            mass = 0
            reading = 0
            while index + 4 < len(vals):
                #print "Mass %d" % mass
                pulse = None
                analog = None
                if vals[index:index+2] != key:
                    raise Exception("missing key")
                index += 2
                # Read the mystery value
                mystery = vals[index+1] * 65536 + vals[index]
                mysteries.append(str(mystery))
                index += 2
                #print "index = %d" % index
                #pprint(vals[index:index+2])
                # Read the pulse count
                #print "Read %d" % reading
                scale = vals[index+1]
                #print "scale = 0x%x" % scale
                pulse = Scale(vals[index], vals[index+1], acf)
                #print "pulse = %s" % pulse
                index += 2
                if vals[index:index+2] != key and not Magic(vals[index+1]):
                    # Read the analog value
                    scale = vals[index+1]
                    analog = Scale(vals[index], vals[index+1], acf)
                    index += 2
                """
                if value is None:
                    raise Exception("scan %d acquisition %d index %d" % (scanNumber, acquisition, index))
                total += value
                n += 1
                """
                if pulse is None and analog is None:
                    pass
                    raise Exception("scan %d acquisition %d index %d" % (scanNumber, acquisition, index))

                if pulse is not None:
                    pulses.append(str(pulse))
                if analog is not None:
                    analogs.append(str(analog))

                reading += 1
                if Magic(vals[index+1]):
                    # End of element
                    #print "end of element"
                    """
                    averages.append(str(total/n))
                    #result.append('X')
                    total = 0.0
                    n = 0
                    """
                    X = vals[index]
                    #result += pulses + analogs + [str(X), '']
                    result += pulses + analogs + ['']
                    mysteries += [''] * (2 + len(analogs))
                    #headers += ["%dp" % elements[mass]] * len(pulses) + ["%da" % elements[mass]] * len(analogs) + ['X', '']
                    headers += ["%dp" % elements[mass]] * len(pulses) + ["%da" % elements[mass]] * len(analogs) + ['']
                    pulses = []
                    analogs = []
                    index += 2
                    mass += 1
                    reading = 0
                # Resynchronize with key. The header/trailer appears to be variable size.
                while index + 4 < len(vals) and vals[index:index+2] != key:
                    #raise Exception("scan %d acquisition %d index %d" % (scanNumber, acquisition, index))
                    #print "index = %d" % index
                    index += 2
                acquisition += 1
            if first == True:
                print >> output, ",".join(headers)
                first = False
            print >> output, ",".join(result)
            #if scanNumber == indexLength:
                #print >> output, "\n" + ",".join(mysteries)
        dat.close()
    output.close()

if __name__ == '__main__':
    main(sys.argv)


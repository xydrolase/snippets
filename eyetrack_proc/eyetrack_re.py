#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """
Processes the data file of a eye tracker, split experiment into trials,
and output summarized eye movement data.
"""

from __future__ import print_function

import re
import sys
import bisect
import numpy

class Trial(list):
    def __init__(self, tid, serial):
        self.time = tid
        self.serial = serial

    def __repr__(self):
        return "{{{0}@{1}}}".format(self.serial, self.time)

class SampleList(list):
    def __init__(self, lst):
        self.is_sorted = False

        super(SampleList, self).__init__(lst)

    def append(self, item):
        self.is_sorted = False

        super(SampleList, self).append(item)

    def extend(self, item):
        self.is_sorted = False

        super(SampleList, self).extend(item)

    def slice(self, start, stop):
        if not self.is_sorted:
            self = SampleList(sorted(self))
            self.is_sorted = True

        lo = bisect.bisect_left(self, start)
        hi = bisect.bisect_right(self, stop)

        #print(start, stop, lo, hi, self[lo:hi])

        return self[lo:hi]

class Event:
    SACCADE = 1
    FIXATION = 2
    types = {'SACC': SACCADE, 'FIX': FIXATION}

    def __init__(self, etype):
        _type = Event.types.get(etype, None)
        if _type is None:
            raise ValueError("Invalid event type.")

        self.event_type = _type
        self.start_time = 0
        self.end_time = 0
        self.blinked = False
        self.assigned_trial = None 

    def __contains__(self, time):
        return self.start_time <= time and self.end_time >= time

    def summarize(self, samples):
        ev_samples = samples.slice(self.start_time, self.end_time)
        self.mean_time = numpy.mean([smp.time for smp in ev_samples])
        self.stats = (
            self.mean_time,
            numpy.mean([smp.x for smp in ev_samples]),
            numpy.std([smp.x for smp in ev_samples]),
            numpy.mean([smp.y for smp in ev_samples]),
            numpy.std([smp.y for smp in ev_samples]),
            numpy.mean([smp.pupil_size for smp in ev_samples]),
            numpy.std([smp.pupil_size for smp in ev_samples])
        )

class Sample:
    __slots__ = ['time', 'x', 'y', 'pupil_size']

    def __init__(self, time, x, y, psize):
        self.time = int(time)
        self.x = float(x)
        self.y = float(y)
        self.pupil_size = float(psize)

    def __repr__(self):
        return "{0} {1} {2} {3}".format(self.time, self.x, self.y, 
                                        self.pupil_size)

    def __eq__(self, t):
        if type(t) is int:
            return self.time == t
        else:
            return self.time == t.time

    def __lt__(self, t):
        if type(t) is int:
            return self.time < t
        else:
            return self.time < t.time

    def __gt__(self, t):
        if type(t) is int:
            return self.time > t
        else:
            return self.time > t.time

    def __le__(self, t):
        if type(t) is int:
            return self.time <= t
        else:
            return self.time <= t.time

    def __ge__(self, t):
        if type(t) is int:
            return self.time >= t
        else:
            return self.time >= t.time

def main():
    regex_segment = re.compile(r"^([SE])([A-Z]+) L\s+(\d+)")
    regex_sample = re.compile(
        r"(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\t\.\.\.")
    regex_trial = re.compile(r"MSG\t(\d+) SYNCTIME")

    trial_serial = 0
    curr_segment = None
    events = []

    lines = None
    with open('s1032_02.asc', 'r') as f:
        lines = [line.rstrip() for line in f]

    print("Processing samples...", file=sys.stderr)

    samples = SampleList([])
    for match in [regex_sample.search(line) for line in lines]:
        if match:
            samples.append(Sample(*match.groups()))

    print("Processing events...", file=sys.stderr)

    # process for fixation and saccade
    for line in lines:
        match = regex_segment.search(line)
        if match:
            start_end, stype, time = match.groups()
            time = int(time)

            if start_end == 'S' and stype in ('SACC', 'FIX'):
                curr_segment = Event(stype)
                curr_segment.start_time = time
            elif start_end == 'E' and stype in ('SACC', 'FIX'):
                if curr_segment and \
                   Event.types.get(stype, None) == curr_segment.event_type:
                    cols = re.split('\s+', line)
                    etime = int(cols[3])
                    curr_segment.end_time = etime
                    curr_segment.summarize(samples)
                    events.append(curr_segment)
                    curr_segment = None
                else:
                    raise TypeError("Inconsistent event segment.")
            elif stype == 'BLINK':
                curr_segment.blinked = True

    print("Processing trials...", file=sys.stderr)

    all_trials = []
    # process for all trails (SYNCTIME)
    for line in lines:
        match = regex_trial.search(line)
        if match:
            trial_time = int(match.group(1))
            trial = Trial(trial_time, trial_serial)

            all_trials.append(trial)

            trial_serial += 1

    # now that all trials, all events and all samples are processed,
    # generate output... 

    print("Generating output...", file=sys.stderr)

    fout = open('eyetracking_re.txt', 'w')
    for ev in events:
        prev_tr = None
        for tr in all_trials:
            if ev.mean_time < tr.time:
                break

            prev_tr = tr

        ev.assigned_trial = prev_tr

        if ev.assigned_trial is None:
            continue

        fout.write("{0} {1} {2}\n".format(
            ev.assigned_trial.serial,
            " ".join([str(st) for st in ev.stats]),
            ev.event_type))

    fout.close()

#def main():
#
#    trials = []
#    saccs = []
#
#    cur_trial = None
#
#    with open('s1032_02.asc', 'r') as f:
#        for line in f:
#            seg_match = regex_segment.search(line)
#            if seg_match:
#                start_end, seg_type = seg_match.groups()
#                if start_end == 'S':
#                    if seg_type == 'SACC': 
#                        saccs.append(SACC(seg_type))
#                    elif seg_type in ('BLINK') and saccs:
#                        saccs[len(saccs)-1].flag |= SACC.FLAG_BLINK
#                else:
#                    if seg_type == 'SACC':
#                        seg = saccs.pop()
#                        if seg.flag == 0 and cur_trial is not None:
#                            cur_trial.extend(seg)
#                        elif seg.flag > 0:
#                            if seg.flag & 0x1:
#                                print(repr(cur_trial), 
#                                        "BLINK between SACC, skip.")
#                            if seg.flag & 0x2:
#                                print(repr(cur_trial), 
#                                        "SYNCTIME between SACC, skip.")
#
#                continue
#
#            trial_match = regex_trial.search(line)
#            if trial_match:
#                if saccs:
#                    saccs[len(saccs)-1].flag |= SACC.FLAG_TRIAL
#                if cur_trial is not None:
#                    trials.append(cur_trial)
#
#                trial_id = int(trial_match.group(1))
#                cur_trial = Trial(trial_id, len(trials))
#
#                continue
#
#            smp_match = regex_sample.search(line)
#            if smp_match:
#                # within SACC segment
#                if saccs:
#                    saccs[len(saccs)-1].append(Sample(*smp_match.groups()))
#                elif cur_trial is not None:
#                    cur_trial.append(Sample(*smp_match.groups()))
#
#        if cur_trial is not None:
#            trials.append(cur_trial)
#
#
#    print("Total trials:", len(trials))
#
#    #print("Trials\n-----\n\n", "\n".join([repr(t) for t in trials]))

if __name__ == "__main__":
    main()

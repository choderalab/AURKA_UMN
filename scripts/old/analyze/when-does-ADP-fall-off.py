import numpy as np
import sys
import math
import os

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
projects = ['11410','11411']
projects = ['11414','11418']
project_dirs = {
    '11410':'%s/../output-1OL5' % local_path,
    '11411':'%s/../output-1OL7' % local_path,
    '11414':'%s/../output-11414' % local_path,
    '11418':'%s/../output-11418' % local_path
}

system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed','11414':'with TPX2'}

if len(sys.argv) == 2:
    this_project = int(sys.argv[1])
    projects = [projects[this_project%len(projects)]]

for i, project in enumerate(projects):
    print('Project %s' % project)
    project_dir = project_dirs[project]
    ADP_bound = np.load('%s/is-ADP-bound.npy' % project_dir)
    for clone, traj in enumerate(ADP_bound):
        if clone%50 == 0:
            print(' RUN%s' % str(clone/50))
            print('  trajectories that begin bound:')
        if traj[0]:
            last_bound_frame = 0
            first_unbound_frame = -1
            found_first_loss = False
            if clone%50 == 14 and clone/50 == 0:
                print('   ** CLONE 14 **')
            for frame_id, frame in enumerate(traj):
                if clone%50 == 14 and clone/50 == 0 and frame_id in range(20):
                    print('    frame%s: %s' % (frame_id, frame))
                if frame:
                    last_bound_frame = frame_id
                if not frame and not found_first_loss:
                    first_unbound_frame = frame_id
                    found_first_loss = True
            if first_unbound_frame - last_bound_frame == 1:
                print('   clone%s loses ADP immediately following frame %s' % (clone%50, last_bound_frame))
            else:
                print('   clone%s loses ADP first in frame %s and finally after frame %s' % (clone%50, first_unbound_frame, last_bound_frame))

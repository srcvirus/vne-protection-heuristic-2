import argparse
import sys
import os
import subprocess


def execute_one_experiment(executable, pn_topology_file, vn_topology_file,
                           location_constraint_file, vnr_root):
    process = subprocess.Popen([executable, '--pn_topology_file=' +
                                pn_topology_file, '--vn_topology_file=' +
                                vn_topology_file,
                                '--location_constraint_file=' + location_constraint_file],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    solution_time = process.stdout.readline().split(':')[1][1:].rstrip('\n')
    with open(os.path.join(vnr_root, "time"), 'w') as f:
        f.write(solution_time)
    print vnr_root + ": " + solution_time


def main():
    # Read the command line arguments and parse them.
    parser = argparse.ArgumentParser(
        description="Script for automating VNE protection experiments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--testcase_root',
        help='Root directory for test cases',
        required=True)
    parser.add_argument(
        '--executable',
        help='Name of the executable file to run',
        required=True)
    parser.add_argument(
        '--input_type',
        help='Type of the input (synthetic/rocketfuel)',
        required=True)
    args = parser.parse_args()
    root = args.testcase_root
    input_type = args.input_type
    executable = './' + args.executable
    subprocess.Popen(['make'], shell=False, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    sn_set = sorted(tuple(os.walk(root))[0][1])
    for sn in sn_set:
        sn_root = os.path.join(root, sn)
        if input_type == "synthetic":
            lnrs = sorted(tuple(os.walk(sn_root))[0][1])
            for lnr in lnrs:
                sn_path = os.path.join(sn_root, lnr)
                sn_topology_file = os.path.join(sn_path, "sn.topo")
                vnr_dir = os.path.join(sn_path, "vnr")
#                print "vnr_dir = " + vnr_dir
                vnr_types = tuple(os.walk(vnr_dir))[0][1]
                for vnr_type in vnr_types:
                    vnr_type_dir = os.path.join(vnr_dir, vnr_type)
                    vnr_instance_dirs = tuple(os.walk(vnr_type_dir))[0][1]
                    for vnr_instance in vnr_instance_dirs:
                        vnr_instance_dir = os.path.join(vnr_type_dir, vnr_instance)
                        vn_topology_file = os.path.join(vnr_instance_dir, "vn.topo")
                        location_constraint_file = os.path.join(vnr_instance_dir, "vnloc")
#                        print sn_topology_file, vn_topology_file, location_constraint_file
                        execute_one_experiment(executable, sn_topology_file, 
                                               vn_topology_file,
                                               location_constraint_file,
                                               vnr_instance_dir)
        elif input_type == "rocketfuel":
            sn_topology_file = os.path.join(sn_root, "sn.topo")
            vnr_dir = os.path.join(sn_root, "vnr")
            vnr_types = tuple(os.walk(vnr_dir))[0][1]
            for vnr_type in vnr_types:
                vnr_type_dir = os.path.join(vnr_dir, vnr_type)
                vnr_instance_dirs = tuple(os.walk(vnr_type_dir))[0][1]
                for vnr_instance in vnr_instance_dirs:
                    vnr_instance_dir = os.path.join(vnr_type_dir, vnr_instance)
                    vn_topology_file = os.path.join(vnr_type_dir, vnr_instance, "vn.topo")
                    location_constraint_file = os.path.join(vnr_type_dir, vnr_instance, "vnloc")
                    print sn_topology_file, vn_topology_file, location_constraint_file
                    execute_one_experiment(executable, sn_topology_file,
                                        vn_topology_file, 
                                        location_constraint_file,
                                        vnr_instance_dir)

if __name__ == "__main__":
    main()

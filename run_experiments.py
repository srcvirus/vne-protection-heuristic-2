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
    args = parser.parse_args()
    root = args.testcase_root
    executable = './' + args.executable
    subprocess.Popen(['make'], shell=False, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    testcases = sorted(tuple(os.walk(root))[0][1])
    for test in testcases:
        path = os.path.join(root, test)
        if os.path.isfile(path + "/sn2.txt"):
            pn_topology_file = os.path.join(path, "sn2.txt")
        else:
            pn_topology_file = os.path.join(path, "sn.txt")
        vnr_dirs = list(os.walk(path))[0][1]
        if len(vnr_dirs) <= 0:
            vn_topology_file = os.path.join(path, "vn.txt")
            location_constraint_file = os.path.join(path, "vnloc.txt")
            execute_one_experiment(executable, pn_topology_file, vn_topology_file,
                                   location_constraint_file, path)
        else:
            for vnr_root in vnr_dirs:
                vnr_types = tuple(os.walk(os.path.join(path, vnr_root)))[0][1]
                for vnr_type in vnr_types:
                    vnr_cases = tuple(os.walk(os.path.join(path, vnr_root,
                                                           vnr_type)))[0][1]
                    for vnr_case in vnr_cases:
                        vnr_r = os.path.join(
                            path,
                            vnr_root,
                            vnr_type,
                            vnr_case)
                        vn_topology_file = os.path.join(vnr_r, "vn.txt")
                        location_constraint_file = os.path.join(
                            vnr_r,
                            "vnloc.txt")
                        execute_one_experiment(executable, pn_topology_file,
                                               vn_topology_file, location_constraint_file, vnr_r)


if __name__ == "__main__":
    main()

import os
import re
from collections import deque
import pytest
import helper


class TestOrthoFinderRuns:
    
    @pytest.mark.order(1)
    def test_orthofinder_core_runs(self, projects, dna_projects, species_tree, of_command_dict, capfd):
        
        core_runs = []
        for category, cmd_dict in of_command_dict.items():
            
            if "assign" in category.lower(): 
                continue 
            
            for name, argstr in cmd_dict.items():
                if not argstr:
                    raise AssertionError(f"Commands for test {name} were not found:\n")
                if "DAN_INPUT" in argstr:
                    argstr = argstr.strip().replace("DNA_INPUT", dna_projects)
                    
                argstr = argstr.strip().replace("INPUT", projects)
                argstr = argstr.strip().replace("SPECIES_TREE", species_tree)
                args = argstr.split()[1:]
                args.extend(["-n", name])
                core_runs.append(args)
                
        for args in core_runs:
            print()
            print(f"TESTING: {name} with command {" ".join(args)}...")
            code, out, err, text = helper._run_main(args, capfd)
            assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

            fatal_patterns = (
                r"(?mi)^\s*ERROR:\s",              
                r"(?m)^\s*Traceback \(most recent call last\):",    
                r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
                r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
            )
            assert not any(re.search(p, text) for p in fatal_patterns), text
            
            
    @pytest.mark.order(2)
    def test_orthofinder_assign_runs(self, projects, assign, species_tree_assign, of_command_dict, capfd):
        
        assign_runs = []
        for category, cmd_dict in of_command_dict.items():
            
            if "core" in category.lower(): 
                continue 
            
            for name, argstr in cmd_dict.items():
                if not argstr:
                    raise AssertionError(f"Commands for test {name} were not found:\n")
                
                argstr = argstr.strip().replace("INPUT", assign)

                results_dir = os.path.join(projects, "OrthoFinder")
                projects_results = helper._find_output_dir(results_dir, "Results_" + name.rsplit("_", 1)[0])
                argstr = argstr.strip().replace("CORE_RESULTS", projects_results)
                argstr = argstr.strip().replace("INPUT", projects)
                argstr = argstr.strip().replace("SPECIES_TREE", species_tree_assign)
                args = argstr.split()[1:]
                args.extend(["-n", name])
                assign_runs.append(args)
                
        for args in assign_runs:
            print()
            print(f"TESTING: {name} with command {" ".join(args)}...")
            code, out, err, text = helper._run_main(args, capfd)
            assert code == 0, f"Exit {code}\n--- stdout ---\n{out}\n--- stderr ---\n{err}"

            fatal_patterns = (
                r"(?mi)^\s*ERROR:\s",              
                r"(?m)^\s*Traceback \(most recent call last\):",    
                r'(?m)^\s*File ".*", line \d+, in \S+',    # Python stack frame lines
                r"(?m)^[A-Za-z_]\w*Error:\s",    # e.g., ValueError:
            )
            assert not any(re.search(p, text) for p in fatal_patterns), text

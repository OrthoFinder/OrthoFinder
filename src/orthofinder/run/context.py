"""The state of one OrthoFinder run, shared by the parts of run.main."""

from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class RunContext:
    log: Any = None                  # utils.logging.Logger of the run (checkpoint.txt)
    current_step: str = "Initialisation"  # what the error handlers of main() say failed
    options: Any = None              # process_args.Options
    prog_caller: Any = None          # program_caller.ProgramCaller
    fasta_dir: Optional[str] = None
    continuation_dir: Optional[str] = None
    results_dir: Optional[str] = None
    state: dict = field(default_factory=dict)  # run_state of the run being continued
    i_og_restart_state: int = 0
    assign_resume: bool = False      # -b on the results of an --assign run

    @contextmanager
    def step(self, name):
        """
        A workflow step: logs "[name] Started", and when it finishes
        "[name] Completed". If it fails, is interrupted or exits, it is not
        completed: current_step stays its name, for the error handlers.
        """
        self.current_step = name
        self.log.step(name, "Started")
        yield
        self.log.step(name, "Completed")
        self.current_step = "Workflow"

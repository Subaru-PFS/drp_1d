import concurrent.futures
import subprocess
import textwrap
import os
import uuid


class Worker:
    def __init__(self, config):
        self.concurrency = config.concurrency
        self.notification_url = config.notification_url

    def run(self, args):
        pass

    def wait_all(self):
        pass


class Local(Worker):
    """Multi-process worker"""
    def __init__(self, config):
        super().__init__(config)
        self.executor = concurrent.futures.ThreadPoolExecutor(self.concurrency)
        self.futures = []

    def run(self, args):
        f = self.executor.submit(subprocess.run, args)
        self.futures.append(f)

    def wait_all(self):
        concurrent.futures.wait(self.futures)


class PBS(Worker):
    """PBS worker"""

    script_template = textwrap.dedent("""\
            #PBS -N {executor_script}
            #PBS -l nodes=1:ppn=1
            #PBS -l walltime=01:00:00
            #PBS -t 1-{jobs}
            cd {workdir}
            {pre_command}
            /usr/bin/env python3 {executor_script} ${{PBS_ARRAYID}} >> {output_dir}/out-{task_id}-${{PBS_ARRAYID}}.txt
            echo "$?" >> {output_dir}/tmp/{task_id}_${{PBS_ARRAYID}}.done
            """)

    def __init__(self, config):
        super().__init__(config)
        self.output_dir = os.path.abspath(config.output_dir)
        self.pre_command = config.pre_command
        self.workdir = os.getcwd()
        self.tasks = []

    def run(self, args):
        self.tasks.append(args)

    def wait_all(self):
        """Submit all jobs to PBS."""

        task_id = uuid.uuid4().hex
        os.makedirs(os.path.join(self.output_dir, 'tmp'), exist_ok=True)
        executor_script = os.path.join(self.output_dir, 'tmp',
                                       'pbs_executor_{}.py'.format(task_id))

        # generate pbs script
        script = PBS.script_template.format(jobs=len(self.tasks),
                                            workdir=self.workdir,
                                            output_dir=self.output_dir,
                                            pre_command=self.pre_command if self.pre_command else '',
                                            executor_script=executor_script,
                                            task_id=task_id)
        pbs_script_name = os.path.join(self.output_dir, 'tmp',
                                       'pbs_script_{}.sh'.format(task_id))
        with open(pbs_script_name, 'w') as pbs_script:
            pbs_script.write(script)

        # generate executor script
        with open(os.path.join(os.path.dirname(__file__),
                               'pbs_executor.py.in'),
                  'r') as f:
            pbs_executor = f.read().format(tasks=self.tasks,
                                           notification_url=self.notification_url)
        with open(executor_script, 'w') as executor:
            executor.write(pbs_executor)

        # run pbs
        result = subprocess.run(['qsub', pbs_script_name])
        assert result.returncode == 0

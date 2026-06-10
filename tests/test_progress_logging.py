import unittest

from fp_tools.utils import utilities


class ProgressLoggingTest(unittest.TestCase):
    def test_progress_milestones_are_coarse_and_include_completion(self):
        seen = []
        for done in range(0, 101):
            pct = utilities.progress_log_percent(done, 100, previous_percent=seen[-1] if seen else None)
            if pct is not None:
                seen.append(pct)

        self.assertEqual(seen, [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

    def test_worker_tqdm_is_disabled_for_normal_verbosity(self):
        self.assertFalse(utilities.show_worker_progress(verbosity=3, total_items=100, is_tty=True))
        self.assertTrue(utilities.show_worker_progress(verbosity=5, total_items=100, is_tty=True))
        self.assertFalse(utilities.show_worker_progress(verbosity=5, total_items=100, is_tty=False))

    def test_monitor_progress_logs_coarse_milestones(self):
        class FakeTask:
            def __init__(self, ready_at, counter):
                self.ready_at = ready_at
                self.counter = counter

            def ready(self):
                return self.counter["value"] >= self.ready_at

        class FakeLogger:
            def __init__(self):
                self.messages = []

            def info(self, message):
                self.messages.append(message)

        counter = {"value": 0}
        tasks = [FakeTask(i, counter) for i in range(1, 101)]
        logger = FakeLogger()
        original_sleep = utilities.time.sleep

        def advance(_seconds):
            counter["value"] += 1

        utilities.time.sleep = advance
        try:
            utilities.monitor_progress(tasks, logger)
        finally:
            utilities.time.sleep = original_sleep

        self.assertEqual(
            logger.messages,
            [f"Progress {pct}%" for pct in range(0, 101, 10)],
        )


if __name__ == "__main__":
    unittest.main()

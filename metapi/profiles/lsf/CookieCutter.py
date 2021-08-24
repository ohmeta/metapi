class CookieCutter:
    """
    Cookie Cutter wrapper
    """

    @staticmethod
    def get_default_threads() -> int:
        return int("1")

    @staticmethod
    def get_default_mem_mb() -> int:
        return int("1024")

    @staticmethod
    def get_log_dir() -> str:
        return "logs/cluster"

    @staticmethod
    def get_default_queue() -> str:
        return ""

    @staticmethod
    def get_lsf_unit_for_limits() -> str:
        return "KB"

    @staticmethod
    def get_unknwn_behaviour() -> str:
        return "wait"

    @staticmethod
    def get_zombi_behaviour() -> str:
        return "ignore"

    @staticmethod
    def get_latency_wait() -> float:
        return float("5")

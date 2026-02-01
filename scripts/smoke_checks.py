def main():
    from metalminer.core import Config
    from metalminer.cli import build_parser

    _ = Config()
    parser = build_parser()
    parser.parse_args([])
    print("Smoke checks passed: imports and CLI parser initialized.")


if __name__ == "__main__":
    main()

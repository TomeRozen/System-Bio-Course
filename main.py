from Game import Game
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Game of Life simulation')

    parser.add_argument('--board_size', '-s', type=int, default=100, help='Size of the board edge')
    parser.add_argument('--begin_alive', '-a', type=int, default=2500, help='Number of cells that begin alive')
    parser.add_argument('--turns', '-t', type=int, default=10000, help='Number of turns to simulate')

    args = parser.parse_args()

    game = Game(args.board_size, args.begin_alive, args.turns)
    game.run()



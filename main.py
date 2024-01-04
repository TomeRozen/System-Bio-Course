from time import sleep
from Board import Board
import numpy as np

if __name__ == '__main__':
    board = Board(100)
    board.random_init(2500)
    board.print_board()


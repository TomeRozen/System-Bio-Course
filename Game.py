import numpy as np

from Board import Board

BOARD_SIZE = 100
BEGIN_ALIVE = 2500
TURNS = 10000


class Game:
    def __init__(self, board_size=BOARD_SIZE, begin_alive=BEGIN_ALIVE, turns=TURNS):
        if begin_alive > board_size**2:
            raise ValueError("Number of cells that begin alive must be smaller than size of board!")
        self.board = Board(board_size)
        self.board.random_init(begin_alive)
        self.turns = turns

    def new_state(self, index: int) -> int:
        """Returns the new state of a cell, using the rules of Game of Life, and the state of its neighbours"""
        # translate index to x and y coordinates of cell
        x, y = np.unravel_index(index, (self.board.get_size(), self.board.get_size()))
        neighbours = self.board.get_num_neighbours(x, y)
        if self.board.get_value(x, y) == 1:
            # if cell is alive
            if neighbours < 2 or neighbours > 3:
                return 0
            else:
                return 1
        else:
            # if cell is dead
            if neighbours == 3:
                return 1
            else:
                return 0


    def turn(self) -> None:

        vec_new_state = np.vectorize(self.new_state)
        # Calculate new state of all cells
        new_board = vec_new_state(np.arange(self.board.get_size()**2))

        self.board.set_board(np.reshape(new_board, (self.board.get_size(), self.board.get_size())))

    def run(self):
        for i in range(self.turns):
            old_board = self.board.get_board().copy()
            self.board.print_board()
            self.turn()
            if self.board.board_sum() == 0:
                break
            if np.array_equal(old_board, self.board.get_board()):
                break

        self.board.print_board()




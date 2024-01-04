import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
class Board:
    def __init__(self, size: int):
        self.size = size
        self.board = np.zeros((size, size), dtype=int)
        self.fig, self.ax = plt.subplots()

    def get_board(self) -> np.ndarray:
        return self.board

    def get_size(self) -> int:
        return self.size

    def set_board(self, board: np.ndarray) -> None:
        self.board = board.copy()

    def set_size(self, size: int) -> None:
        self.size = size

    def get_value(self, x: int, y: int) -> int:
        return self.board[x, y]

    def set_value(self, x: int, y: int, value: int) -> None:
        self.board[x, y] = value

    def random_init(self, num_on: int) -> None:
        self.reset_board()
        on_indices = np.random.choice(self.size**2, num_on, replace=False)
        self.board[np.unravel_index(on_indices, (self.size, self.size))] = 1


    def get_num_neighbours(self, x: int, y: int) -> int:
        """Returns the number of neighbours of a cell"""
        neighbours = 0
        if x > 0:
            neighbours += self.board[x-1, y]
        if x < self.size - 1:
            neighbours += self.board[x+1, y]
        if y > 0:
            neighbours += self.board[x, y-1]
        if y < self.size - 1:
            neighbours += self.board[x, y+1]
        return neighbours

    def reset_board(self) -> None:
        self.set_board(np.zeros((self.size, self.size), dtype=int))

    def print_board(self) -> None:
        """prints the board nicely to the console"""
        # plt.clf()
        self.ax.imshow(self.board, cmap=plt.cm.binary)
        plt.pause(10)

    def board_sum(self):
        return np.sum(self.board)




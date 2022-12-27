import src.fat_tait as ft
import sys

if __name__ == '__main__':
    gauss_4_105 = sys.argv[1]
    # Example: 'O1-U2-O3-U1-O4+U3-O2-U4+'
    g, h = ft.get_fat_tait_graphs(gauss_4_105)

    print('First fat Tait graph:')
    print('[edge list]')
    print(g.edge_list)
    print('[cyclic ordering dictionary]')
    print(g.cord_dict)
    print()
    print('Second fat Tait graph:')
    print('[edge list]')
    print(h.edge_list)
    print('[cyclic ordering dictionary]')
    print(h.cord_dict)

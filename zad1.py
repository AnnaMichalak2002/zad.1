import sys
import numpy as np    #do macierzy kontaktów
from Bio import PDB     #do analizy struktury białka z pliku PDB.
import matplotlib.pyplot as plt  #by wykres mieć

def calculate_contact_map(structure, threshold=8.0): #próg odl = 8
    atom_list = list(structure.get_atoms()) # konwertujemy metodę w bibliotece Biopython, która zwraca generator zawierający wszystkie atomy w strukturze białka na liste - mamy liste atomów
    num_atoms = len(atom_list) # ilość atomów
    
    contact_map = np.zeros((num_atoms, num_atoms), dtype=int) #inicjujemy macierz kontaktów wypelnioną zerami z biblioteki numpy o wartosciach calkowitoliczbowych
    
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):  # iterujemy przez każdą parę atomów CA w białku
            distance = atom_list[i] - atom_list[j]  #obliczamy dystans
            if distance < threshold:         #jesli dystans mniejszy niz threshold to dajemy 1 w mapie kontaktów
                contact_map[i, j] = 1
                contact_map[j, i] = 1
    
    return contact_map

def plot_contact_map(contact_map):    #tu wykorzystamy matplotliba do stworzenia wizualizacji mapy kontaktów za pomocą wykresu typu heatmap
    plt.imshow(contact_map, cmap='viridis', interpolation='none') #funkcji imshow z biblioteki Matplotlib do wyświetlenia macierzy kontaktów jako obrazu,  kolory z palety virdis, brak interpolacji
    plt.title('Contact Map') #tytuł
    plt.xlabel('Residue Index') # oś X odnosi się do indeksów reszt białka w macierzy kontaktów.
    plt.ylabel('Residue Index')  # oś Y  odnosi się do indeksów reszt białka w macierzy kontaktów.
    plt.colorbar()  # Dodaje pasek kolorów, który przyporządkowuje kolorystykę z wykresu do odpowiadających im wartości w macierzy kontaktów.
    plt.show()   #wyświetla

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]

    # Wczytujemy strukturę białka
    pdb_parser = PDB.PDBParser(QUIET=True)  #tworzymy obiekt parsera plików PDB z modułu Biopython Bio.PDB Argument QUIET=True wskazuje, że parser działa w trybie cichym czyli nie wyświetla ostrzeżeń ani komunikatów odczytu podczas przetwarzania pliku.
    structure = pdb_parser.get_structure('protein', pdb_file)  #wykorzystujemy parser PDB do odczytania struktury białka z pliku PDB.

    # Wyznaczamy mapę kontaktów
    contact_map = calculate_contact_map(structure)

    # Wyświetlamy mapę kontaktów
    plot_contact_map(contact_map)

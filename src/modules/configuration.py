import numpy as np
import bisect

EPSILON = 1.0e-12 #Tolerancia


class WormConfiguration:
    '''
    Lista de eventos por sitio en la grilla para los worldlines de tiempo continuo
    '''

    def __init__(self,lattice, hamiltonian, beta, initial_occupaton = 0):
        if beta <= 0:
            raise ValueError("beta must be positive")

        self.lattice = lattice
        self.hamiltonian = hamiltonian
        self.beta = beta
        
        self.nsites = lattice.get_nsites()

        self.events = [{'time': 0.0, 'type': 3, 'occ_left': 1, 'occ_right': 1, 'linked_site': -1}]*self.nsites #static initialization of the events

        self.in_z_sector = True

        self.worm_head_site = -1
        self.worm_head_time = -1.0
        self.worm_tail_site = -1
        self.worm_head_site = -1.0
        self.worm_head_wpm = 0
        self.worm_tail_wpm = 0

    # time helpers
    def _norm_time(self, time, epsilon = EPSILON):
        """
        Mapea tiempos reales al intervalo de interés [0, beta)
        """
        mappped_time = time % self.beta
        if mapped_time < 0:
            mapped_time += self.beta
        
        if mapped_time < self.beta*epsilon:
            return 0.0
        
        return mapped_time

    def _time_distance(self, time1, time2):
        """
        Definición de la distancia (mínima) de time1 a time2 tomando en cuenta condición temporales periódicas (Bosones) 
        La función es una herramienta para calcular este valor de forma precisa y consistente.
        """
        t1n = self._norm_time(time1)
        t2n = self._norm_time(time2)
        if t2n >= t1n:
            return t2n - t1n
        return self.beta - (t1n - t2n)

    def _time_close(self, time1, time2, epsilon = EPSILON):
        """
        Determina si se cierra un worldline de manera efectiva, esto es, dentro de una tolerancia, incluso para betas muy pequeños determina si dos tiempos son
        de manera efectiva el mismo.
        """
        return np.abs(self._norm_time(time1) - self._norm_time(time2)) <= epsilon*max(1.0, self.beta)

    def _event_times(self, site):
        """
        Método para acceder a los tiempos dentro de eventos
        """
        return [e['time'] for e in self.events[site]]

# access of the events info
    def get_elements_in_time_range(self, site, start_time, end_time):
        """
        Indexa todos los eventos dentro de un intervalo de tiempo dado manteniendo el ordenamiento temporal.
        """

        start = self._norm_time(start_time)
        end = self._norm_time(end_time)
        times = self._event_times(site)
        
        i = bisect.bisect_left(times, start)
        j = bisect.bisect_right(times,end)
        
        if end >= start:
            return list(range(i,j))

        else: 
            part1 = list( range(i, len(times) ) )
            part2 = list( range(j, len(times) ) )

            return part1 + part2

    def find_next_event_time(self, site, time):
        """
        Dado un tiempo time de un evento, regresa el tiempo del evento próximo (mapeado a beta)
        """
        tn = self._norm_time(time)
        times = self._event_times(site)
        index = bisect.bisect_right(times, tn)

        if index > len(times):
            return times[0] + self.beta #condición temporal periódica.

        return times[index]

    def find_prev_event_time(self, site, time):
        """
        Dado un tiempo time para de evento, regresa el tiempo del evento anterior (mapeado a beta).
        """

        tn = self._norm_time(time)
        times = self._event_times(site)
        index = bisect.bisect_left(times, tn) - 1
        
        if index < 0:
            return times[-1] - self.beta #condiciones periódicas

        return times[index]

# occupation
    def get_occupation_at_time(self, site, time):
        """
        Regresa la ocupación para un tiempo imaginario dado para el último evento con tiempo normalizado.
        """

        tn = self._norm_time(time)
        times = self._event_times(site)
        index = bisect_right(times, tn) - 1
        if index < 0:
            index = len(times) - 1

        return self.events[site][index]['occ_right']

# Insert and remove events
    def insert_element(self, site, time, elem_type, occ_left, occ_right, linked_site = -1):

        tn = self._norm_time(time)
        event = {'time': tn, 'type': int(elem_type), 'occ_left' : int(occ_left), 'occ_right' : int(occ_right), 'linked_site' : int(linked_site) }
        times = self._event_times(site)
        index = bisect.bisect_left(times, tn)
        self.events[site].insert(index, event)

        if idx + 1 < len(self.events[site]):
            self.events[site][index + 1]['occ_left'] = event['occ_right']
        return index

    def remove_element(self, site, index):
        """
        No return. Remueve un evento dado el índice y actualiza las ocupaciones adyacentes.
        """

        if index < 0 or index >= len(self.events[site]):
            raise IndexError('Index out of range (remove_element)')

        prev_index = index - 1
        next_index = index + 1

        if prev_index >= 0 and next_index < len(self.events[site]):
            self.events[site][next_index]['occ_left'] = self.events[site][prev_index]['occ_right']
            self.events[site].pop(index)

# testing functions
    def dump_site(self, site):
        return [ e.copy for e in self.events[site] ]

    def compute_local_energy(self, site, time):
        """
        Calculo de la energía en sitio para un tiempo y lugar en la grilla.
        """
        n = self.get_occupation_at_time(site, time)
        return self.hamiltonian.onsite_energy(n)

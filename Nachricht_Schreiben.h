/*
 * Nachricht_Schreiben.h
 *
 *  Created on: 20.05.2010
 *      Author: martin
 */

#ifndef NACHRICHT_SCHREIBEN_HH_
#define NACHRICHT_SCHREIBEN_HH_

inline void Nachricht_Schreiben(string Nachricht, int Priority, int Prioritylevel)
{
	/*
	 * Diese Funktion soll das implementieren von Nachrichten in einem Programm
	 * vereinfachen.  Damit ein Programm schnell läuft, dürfen nicht unentwegt
	 * Nachrichten geschrieben werden.
	 * Um ein Programm aber
	 * a) zu debuggen
	 * b) Warnungen auswerfen zu lassen, falls man vermutet, dass es falsch
	 *    verwendet wird(z.B. Laden einer garnicht vorhandenen Datei etc.) muss
	 *    man hin und wieder Text ausschreiben.
	 * Die Nachrichten aus a) sollen später garnicht mehr gezeigt werden.
	 * Bei den Nachrichten für b) kann man noch Abstufungen treffen, welche
	 * dieser Nachrichten gezeigt werden.
	 * Zur Nomenklatur:
	 * Prioritylevel ist das Mindestlevel der Anzuzeigenden Nachrichten...
	 * hier ist es Sinnvoll dieses als globale Variable
	 * im Hauptprogramm festzulegen
	 * Priority ist die individulle Wichtigkeit der Nachricht auf einer
	 * frei wählbaren Skala...z.B. 0 garnicht wichtig bis 10 Immer zeigen
	 * Nachricht ...ist selbsterklärend
	 * Lange Rede, kurzer Sinn los gehts....
	  */
	if (Priority >= Prioritylevel)
		cerr << Nachricht << "\n";
	// Ja das wars schon....Extrem lange Funktion, nicht?
	// Jetzt kann man noch eine Richtlinie festlegen Welcher Nachricht,
	// welche Priorität haben soll
	/*Priorität Nachricht
	  10         Alles was zum Programmabbruch führt
	   9          Hinweise auf falsche Datenfiles
	   3          Allgemeine durchs Programm führende Hinweise
	   2          Debugging Hinweis, aktuelles Problem
	   1          Debugging Hinweis, nicht aktuelles Problem
	*/
}


#endif /* NACHRICHT_SCHREIBEN_HH_ */

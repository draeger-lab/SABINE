/*
    SABINE predicts binding specificities of transcription factors.
    Copyright (C) 2009 ZBIT, University of Tübingen, Johannes Eichner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package gui;
import java.io.PrintStream;
import java.util.Calendar;

public class MessageProcessor implements MessageListener {

	private PrintStream err;
	private PrintStream out;

	/**
	 * 
	 */
	public MessageProcessor() {
		this(System.out, System.err);
	}

	/**
	 * 
	 * @param out
	 * @param err
	 */
	public MessageProcessor(PrintStream out, PrintStream err) {
		this.out = out;
		this.err = err;
	}

	public PrintStream getErrorStream() {
		return err;
	}

	public PrintStream getOutStream() {
		return out;
	}

	public static final String getTime() {
		Calendar c = Calendar.getInstance();
		StringBuffer sb = new StringBuffer();
		sb.append(twoDigits(c.get(Calendar.YEAR)));
		sb.append('-');
		sb.append(twoDigits(c.get(Calendar.MONTH) + 1));
		sb.append('-');
		sb.append(twoDigits(c.get(Calendar.DATE)));
		sb.append("T");
		sb.append(twoDigits(c.get(Calendar.HOUR_OF_DAY)));
		sb.append(':');
		sb.append(twoDigits(c.get(Calendar.MINUTE)));
		sb.append(':');
		sb.append(twoDigits(c.get(Calendar.SECOND)));
		return sb.toString();
	}

	private static final String twoDigits(int digit) {
		if (digit < 10) {
			StringBuffer sb = new StringBuffer();
			sb.append(Character.valueOf('0'));
			sb.append(Integer.toString(digit));
			return sb.toString();
		}
		return Integer.toString(digit);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see quick.io.MessageListener#logError(java.lang.Object)
	 */
	public void logError(Object message) {
		//err.print('[');
		//err.print(getTime());
		//err.print("]\t");
		err.println(message);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see quick.io.MessageListener#logMessage(java.lang.Object)
	 */
	public void logMessage(Object message) {
		out.print('[');
		out.print(getTime());
		out.print("]\t");
		out.println(message);
	}

	public void setErrorStream(PrintStream err) {
		this.err = err;
	}

	public void setOutStream(PrintStream out) {
		this.out = out;
	}
}


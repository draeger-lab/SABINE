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

package help;

import java.text.DecimalFormat;
import java.util.Calendar;

public class TimeStampGenerator {

	public String getTimeStamp() {
		
		DecimalFormat fmt = new DecimalFormat();
		fmt.setMaximumIntegerDigits(2);
		fmt.setMinimumIntegerDigits(2);
		
		// get current time and date
		Calendar cal = Calendar.getInstance ();
		String curr_date = (fmt.format(cal.get(Calendar.DAY_OF_MONTH)) + "." + 
							fmt.format((cal.get(Calendar.MONTH) + 1)) + "." + 
							cal.get(Calendar.YEAR));
		String curr_time = (fmt.format(cal.get(Calendar.HOUR_OF_DAY)) + ":" +
							fmt.format(cal.get(Calendar.MINUTE)));
		
		String time_stamp = curr_date + "_" + curr_time;
		
		return time_stamp;
	}
}


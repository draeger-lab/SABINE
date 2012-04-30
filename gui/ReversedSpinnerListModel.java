package gui;

import javax.swing.SpinnerListModel;
import javax.swing.SpinnerModel;

public class ReversedSpinnerListModel extends SpinnerListModel {

	private static final long serialVersionUID = 1L;
	
	Object firstValue, lastValue;
    SpinnerModel linkedModel = null;

    public ReversedSpinnerListModel(Object[] values) {
        super(values);
        firstValue = values[0];
        lastValue = values[values.length - 1];
    }

    public void setLinkedModel(SpinnerModel linkedModel) {
        this.linkedModel = linkedModel;
    }

    public Object getNextValue() {
        Object value = super.getPreviousValue();
        return value;
    }

    public Object getPreviousValue() {
        Object value = super.getNextValue();
        return value;
    }
}


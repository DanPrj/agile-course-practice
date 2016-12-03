package ru.unn.agile.personalfinance.view.controllers;

import javafx.beans.binding.Bindings;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Button;
import javafx.scene.control.TextField;
import javafx.util.converter.CurrencyStringConverter;
import ru.unn.agile.PersonalFinance.ViewModel.AccountViewModel;
import ru.unn.agile.personalfinance.view.ViewModelService;
import ru.unn.agile.personalfinance.view.WindowsManager;

import java.net.URL;
import java.util.ResourceBundle;

public class AddAccountController implements Initializable {
    private final AccountViewModel account =
            new AccountViewModel(ViewModelService.getViewModel());

    @FXML
    private TextField nameField;

    @FXML
    private TextField balanceField;

    @FXML
    private Button addButton;

    public void handleAddButton(final ActionEvent actionEvent) {
        account.save();
        WindowsManager.getInstance().goBack();
    }

    public void handleCancelButton(final ActionEvent actionEvent) {
        WindowsManager.getInstance().goBack();
    }

    @Override
    public void initialize(final URL location, final ResourceBundle resources) {
        setUpBindings();
    }

    private void setUpBindings() {
        /* nameField.text <-> account.name */
        Bindings.bindBidirectional(nameField.textProperty(), account.nameProperty());

        /* balanceField.text <-> account.balance */
        Bindings.bindBidirectional(
                balanceField.textProperty(),
                account.balanceProperty(),
                new CurrencyStringConverter());

        /* account.isAbleToSave -> addButton.disabled */
        addButton.disableProperty().bind(account.isAbleToSaveProperty().not());
    }
}

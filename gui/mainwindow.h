#pragma once

#include <QMainWindow>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QPushButton>
#include <QLineEdit>
#include <QCheckBox>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QTextEdit>
#include <QProcess>
#include <QComboBox>
#include <QRadioButton>
#include <QButtonGroup>
#include <QStackedWidget>
#include <QSettings>
#include <QProgressBar>
#include <QSpinBox>
#include <QScrollArea>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override = default;

private slots:
    void onRunClicked();
    void onBrowseScript();
    void onProcessOutput();
    void onProcessFinished(int exitCode);

private:
    void setupUI();
    void buildAndWriteConfig();
    QString scriptPath() const;

    // Script path
    QLineEdit *m_scriptPathEdit;

    // --- Simulation parameters ---
    QLineEdit *m_proteinNameEdit;
    QDoubleSpinBox *m_temperatureSpin; // K
    QDoubleSpinBox *m_simTimeSpin;     // ns
    QComboBox *m_debugLevelCombo;      // int
    //QCheckBox *m_hisManualCheck;       // yes/no
    QComboBox *m_waterModelCombo;      // What water model is used?
    QSpinBox *m_nrunsSpin;             // How many simulations to run?
    QCheckBox *m_disulfideCheck;       // pdb2gmx -ss
    QComboBox *m_forceFieldCombo;      // force field

    QCheckBox  *m_boxManualCheck;
    QWidget    *m_boxManualWidget;   // container shown/hidden
    QLineEdit  *m_boxDimEdit;
    QComboBox  *m_cellShapeCombo;

    QRadioButton *m_ionModeSalt;
    QRadioButton *m_ionModeManual;
    QButtonGroup *m_ionModeGroup;
    QStackedWidget *m_ionStack;     // page 0 = salt, page 1 = manual
    QDoubleSpinBox *m_saltConcSpin; // already exists, move into stack
    QSpinBox *m_posIonsSpin;
    QSpinBox *m_negIonsSpin;
    QCheckBox *m_magnesiumCheck;

    // Run control
    QPushButton *m_runButton;
    QPushButton *m_browseButton;

    // Progress
    QProgressBar *m_progressBar;    
    int m_totalSteps = 0;

    // Output log
    QTextEdit *m_logOutput;

    // Process
    QProcess *m_process;
};

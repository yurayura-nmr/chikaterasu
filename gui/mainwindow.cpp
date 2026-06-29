#include "mainwindow.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include <QTextStream>
#include <QSizePolicy>
#include <QFont>
#include <QFrame>
#include <QScrollBar>
#include <QComboBox>
#include <QRegularExpression>

// ── colour / style tokens ──────────────────────────────────────────────────
static const char *kStyleSheet = R"(
QMainWindow, QWidget#central {
    background: #0f1117;
}
QGroupBox {
    color: #a0aec0;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    border: 1px solid #2d3748;
    border-radius: 6px;
    margin-top: 14px;
    padding-top: 8px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    top: 2px;
}
QLabel {
    color: #cbd5e0;
    font-size: 13px;
}
QLabel#sectionTitle {
    color: #e2e8f0;
    font-size: 22px;
    font-weight: 700;
}
QLabel#subtitle {
    color: #4a5568;
    font-size: 12px;
}
QLineEdit, QDoubleSpinBox, QSpinBox {
    background: #1a202c;
    border: 1px solid #2d3748;
    border-radius: 5px;
    color: #e2e8f0;
    font-size: 13px;
    padding: 5px 8px;
    min-height: 28px;
}
QLineEdit:focus, QDoubleSpinBox:focus {
    border-color: #4299e1;
}
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
    background: #2d3748;
    border: none;
    width: 18px;
}
QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {
    background: #4a5568;
}
QPushButton#runButton {
    background: #3182ce;
    color: #ffffff;
    border: none;
    border-radius: 6px;
    font-size: 14px;
    font-weight: 600;
    padding: 10px 28px;
    min-height: 38px;
}
QPushButton#runButton:hover  { background: #4299e1; }
QPushButton#runButton:pressed{ background: #2b6cb0; }
QPushButton#runButton:disabled{ background: #2d3748; color: #4a5568; }
QPushButton#browseButton {
    background: #2d3748;
    color: #a0aec0;
    border: 1px solid #4a5568;
    border-radius: 5px;
    font-size: 12px;
    padding: 5px 12px;
    min-height: 28px;
}
QPushButton#browseButton:hover { background: #4a5568; color: #e2e8f0; }
QTextEdit {
    background: #0d1117;
    border: 1px solid #2d3748;
    border-radius: 5px;
    color: #68d391;
    font-family: "Menlo", "Courier New", monospace;
    font-size: 12px;
    padding: 6px;
}
QFrame#divider {
    color: #2d3748;
}
QProgressBar {
    background: #1a202c;
    border: 1px solid #2d3748;
    border-radius: 5px;
    color: #e2e8f0;
    text-align: center;
    font-size: 12px;
    min-height: 22px;
}
QProgressBar::chunk {
    background: #3182ce;
    border-radius: 4px;
}
    QScrollArea {
    background: #0f1117;
}
QScrollBar:vertical {
    background: #1a202c;
    width: 8px;
    border-radius: 4px;
}
QScrollBar::handle:vertical {
    background: #4a5568;
    border-radius: 4px;
    min-height: 20px;
}
QScrollBar::handle:vertical:hover {
    background: #718096;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}
QPushButton#stopButton {
    background: #742a2a;
    color: #fed7d7;
    border: none;
    border-radius: 6px;
    font-size: 14px;
    font-weight: 600;
    padding: 10px 20px;
    min-height: 38px;
}
QPushButton#stopButton:hover   { background: #9b2c2c; }
QPushButton#stopButton:pressed { background: #63171b; }
QPushButton#stopButton:disabled{ background: #2d3748; color: #4a5568; }
)";

// ── constructor ────────────────────────────────────────────────────────────
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), m_process(new QProcess(this))
{
    setWindowTitle("Chikaterasu  —  GROMACS launcher");
    setMinimumWidth(560);
    setStyleSheet(kStyleSheet);
    setupUI();
    resize(620, 860);

    connect(m_process, &QProcess::readyReadStandardOutput,
            this, &MainWindow::onProcessOutput);
    connect(m_process, &QProcess::readyReadStandardError,
            this, &MainWindow::onProcessOutput);
    connect(m_process, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
            this, [this](int exitCode, QProcess::ExitStatus)
            { onProcessFinished(exitCode); });
}

// ── UI construction ────────────────────────────────────────────────────────
void MainWindow::setupUI()
{
    // Outer widget holds the scroll area
    auto *outerWidget = new QWidget(this);
    auto *outerLayout = new QVBoxLayout(outerWidget);
    outerLayout->setContentsMargins(0, 0, 0, 0);
    setCentralWidget(outerWidget);

    // Scroll area
    auto *scroll = new QScrollArea(this);
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    outerLayout->addWidget(scroll);

    // Inner widget is what actually holds all the content
    auto *central = new QWidget();
    central->setObjectName("central");
    scroll->setWidget(central);

    auto *root = new QVBoxLayout(central);
    root->setContentsMargins(24, 20, 24, 20);
    root->setSpacing(16);

    // ── Header ──────────────────────────────────────────────────────────
    auto *titleLabel = new QLabel("Chikaterasu", this);
    titleLabel->setObjectName("sectionTitle");

    auto *subtitleLabel = new QLabel("GROMACS automated MD launcher", this);
    subtitleLabel->setObjectName("subtitle");

    root->addWidget(titleLabel);
    root->addWidget(subtitleLabel);

    // thin divider
    auto *div = new QFrame(this);
    div->setObjectName("divider");
    div->setFrameShape(QFrame::HLine);
    root->addWidget(div);

    // ── Script path ─────────────────────────────────────────────────────
    auto *scriptGroup = new QGroupBox("Script", this);
    auto *scriptLayout = new QHBoxLayout(scriptGroup);
    scriptLayout->setContentsMargins(12, 16, 12, 12);

    m_scriptPathEdit = new QLineEdit(this);
    m_scriptPathEdit->setPlaceholderText("/path/to/chikaterasu.sh");

    m_browseButton = new QPushButton("Browse…", this);
    m_browseButton->setObjectName("browseButton");
    connect(m_browseButton, &QPushButton::clicked, this, &MainWindow::onBrowseScript);

    scriptLayout->addWidget(m_scriptPathEdit, 1);
    scriptLayout->addWidget(m_browseButton);
    root->addWidget(scriptGroup);

    // ── Simulation parameters ────────────────────────────────────────────
    auto *simGroup = new QGroupBox("Simulation parameters", this);
    auto *form = new QFormLayout(simGroup);
    form->setContentsMargins(12, 16, 12, 12);
    form->setSpacing(10);
    form->setLabelAlignment(Qt::AlignRight | Qt::AlignVCenter);

    // Protein / molecule name
    m_proteinNameEdit = new QLineEdit(this);
    m_proteinNameEdit->setText("1UBQ");
    m_proteinNameEdit->setPlaceholderText("e.g. 1UBQ or ATP");
    form->addRow("Molecule / PDB name:", m_proteinNameEdit);

    m_nrunsSpin = new QSpinBox(this);
    m_nrunsSpin->setRange(1, 20);
    m_nrunsSpin->setValue(1);
    form->addRow("Number of runs:", m_nrunsSpin);

    // Temperature
    m_temperatureSpin = new QDoubleSpinBox(this);
    m_temperatureSpin->setRange(0.0, 600.0);
    m_temperatureSpin->setSingleStep(10.0);
    m_temperatureSpin->setDecimals(1);
    m_temperatureSpin->setValue(300.0);
    m_temperatureSpin->setSuffix("  K");
    form->addRow("Temperature:", m_temperatureSpin);

    // Select Force field from Combo box
    m_forceFieldCombo = new QComboBox(this);
    m_forceFieldCombo->addItem("AMBER99SB-ILDN (recommended)", "amber99sb-ildn");
    m_forceFieldCombo->addItem("AMBER99SB", "amber99sb");
    m_forceFieldCombo->addItem("AMBER99", "amber99");
    m_forceFieldCombo->addItem("AMBER96", "amber96");
    m_forceFieldCombo->addItem("AMBER94", "amber94");
    m_forceFieldCombo->addItem("AMBER03", "amber03");
    m_forceFieldCombo->addItem("AMBER14SB", "amber14sb");
    m_forceFieldCombo->addItem("AMBER19SB", "amber19sb");
    m_forceFieldCombo->addItem("AmberGS", "amberGS");
    m_forceFieldCombo->addItem("CHARMM27", "charmm27");
    m_forceFieldCombo->addItem("GROMOS43A1", "gromos43a1");
    m_forceFieldCombo->addItem("GROMOS43A2", "gromos43a2");
    m_forceFieldCombo->addItem("GROMOS45A3", "gromos45a3");
    m_forceFieldCombo->addItem("GROMOS53A5", "gromos53a5");
    m_forceFieldCombo->addItem("GROMOS53A6", "gromos53a6");
    m_forceFieldCombo->addItem("GROMOS54A7", "gromos54a7");
    m_forceFieldCombo->addItem("OPLS-AA", "oplsaa");
    m_forceFieldCombo->setCurrentIndex(0); // AMBER99SB-ILDN as default
    form->addRow("Force field:", m_forceFieldCombo);

    // Select water model from Combo box
    m_waterModelCombo = new QComboBox(this);
    m_waterModelCombo->addItem("TIP3P", "tip3p");
    m_waterModelCombo->addItem("TIP4P", "tip4p");
    m_waterModelCombo->addItem("TIP5P", "tip5p");
    m_waterModelCombo->addItem("SPC", "spc");
    m_waterModelCombo->addItem("SPC/E", "spce");
    m_waterModelCombo->addItem("None", "none");
    m_waterModelCombo->setCurrentIndex(1); // TIP4P default but do not use it for shear-flow
    form->addRow("Water model:", m_waterModelCombo);

    // Debug level
    m_debugLevelCombo = new QComboBox(this);
    m_debugLevelCombo->addItem("[MD] Full production run (NPT)", 0);
    m_debugLevelCombo->addItem("[1] Topology generation", 1);
    m_debugLevelCombo->addItem("[2] Solvation", 2);
    m_debugLevelCombo->addItem("[3] Counterions restraints", 3);
    m_debugLevelCombo->addItem("[4] Energy minimization", 4);
    m_debugLevelCombo->addItem("[5] NVT equilibration", 5);
    m_debugLevelCombo->addItem("[6] NPT equilibration", 6);
    m_debugLevelCombo->setCurrentIndex(0); // default: full run
    form->addRow("Stop after:", m_debugLevelCombo);

    // Histidine - Let gromacs handle it or specify?
    // m_hisManualCheck = new QCheckBox("Manually specify histidine protonation", this);
    // m_hisManualCheck->setChecked(false);
    // form->addRow("Histidine:", m_hisManualCheck);

    m_disulfideCheck = new QCheckBox("Auto-detect disulfide bridges  (-ss)", this);
    m_disulfideCheck->setChecked(false);
    m_disulfideCheck->setToolTip("Runs pdb2gmx with -ss for automatic SS detection.\n"
                                 "For interactive chain/merge selection, use the CLI directly.");
    form->addRow("Disulfide bonds:", m_disulfideCheck);

    // === Box configuration ===
    auto *boxGroup = new QGroupBox("Box configuration", this);
    auto *boxVbox = new QVBoxLayout(boxGroup);
    boxVbox->setContentsMargins(12, 16, 12, 12);
    boxVbox->setSpacing(8);

    m_boxManualCheck = new QCheckBox("Manually specify box dimensions", this);
    m_boxManualCheck->setChecked(false);
    boxVbox->addWidget(m_boxManualCheck);

    // === Shear flow (Rheo-MD) ===
    auto *shearGroup = new QGroupBox("Shear flow (Rheo-MD)", this);
    auto *shearVbox = new QVBoxLayout(shearGroup);
    shearVbox->setContentsMargins(12, 16, 12, 12);
    shearVbox->setSpacing(8);

    m_shearCheck = new QCheckBox("Enable shear flow", this);
    m_shearCheck->setChecked(false);
    m_shearCheck->setToolTip("Disables pressure coupling and applies a deforming box.\n"
                             "SPC/E water is recommended for stability under shear.");
    shearVbox->addWidget(m_shearCheck);

    auto *shearRateWidget = new QWidget(this);
    auto *shearForm = new QFormLayout(shearRateWidget);
    shearForm->setContentsMargins(0, 4, 0, 0);

    m_shearRateSpin = new QDoubleSpinBox(this);
    m_shearRateSpin->setRange(0.0, 1.0);
    m_shearRateSpin->setSingleStep(0.001);
    m_shearRateSpin->setDecimals(4);
    m_shearRateSpin->setValue(0.01);
    m_shearRateSpin->setSuffix("  nm/ps");
    shearForm->addRow("Shear rate:", m_shearRateSpin);

    shearRateWidget->setVisible(false);
    shearVbox->addWidget(shearRateWidget);
    root->addWidget(shearGroup);

    connect(m_shearCheck, &QCheckBox::toggled, shearRateWidget, &QWidget::setVisible);

    // Container widget — hidden by default
    m_boxManualWidget = new QWidget(this);
    auto *boxForm = new QFormLayout(m_boxManualWidget);
    boxForm->setContentsMargins(0, 4, 0, 0);
    boxForm->setSpacing(8);

    m_boxDimEdit = new QLineEdit(this);
    m_boxDimEdit->setText("8.00   8.00   8.00");
    m_boxDimEdit->setPlaceholderText("x   y   z  (nm)");
    boxForm->addRow("Box dimensions (nm):", m_boxDimEdit);

    m_cellShapeCombo = new QComboBox(this);
    m_cellShapeCombo->addItem("Triclinic", "triclinic");
    m_cellShapeCombo->addItem("Cubic", "cubic");
    m_cellShapeCombo->addItem("Dodecahedron", "dodecahedron");
    m_cellShapeCombo->addItem("Octahedron", "octahedron");
    m_cellShapeCombo->setCurrentIndex(0);
    boxForm->addRow("Cell shape:", m_cellShapeCombo);

    m_boxManualWidget->setVisible(false);
    boxVbox->addWidget(m_boxManualWidget);

    root->addWidget(boxGroup);

    connect(m_boxManualCheck, &QCheckBox::toggled,
            m_boxManualWidget, &QWidget::setVisible);

    // === Ion Configuration ===
    auto *ionGroup = new QGroupBox("Ion configuration", this);
    auto *ionVbox = new QVBoxLayout(ionGroup);
    ionVbox->setContentsMargins(12, 16, 12, 12);
    ionVbox->setSpacing(8);

    // Radio buttons
    auto *radioRow = new QHBoxLayout;
    m_ionModeSalt = new QRadioButton("Salt concentration", this);
    m_ionModeManual = new QRadioButton("Manual ion count", this);
    m_ionModeSalt->setChecked(true);
    m_ionModeGroup = new QButtonGroup(this);
    m_ionModeGroup->addButton(m_ionModeSalt, 0);
    m_ionModeGroup->addButton(m_ionModeManual, 1);
    radioRow->addWidget(m_ionModeSalt);
    radioRow->addWidget(m_ionModeManual);
    radioRow->addStretch();
    ionVbox->addLayout(radioRow);

    // Stacked widget — page 0: concentration, page 1: manual counts
    m_ionStack = new QStackedWidget(this);

    // Page 0
    auto *saltPage = new QWidget;
    auto *saltForm = new QFormLayout(saltPage);
    saltForm->setContentsMargins(0, 4, 0, 0);
    m_saltConcSpin = new QDoubleSpinBox;
    m_saltConcSpin->setRange(0.0, 2.0);
    m_saltConcSpin->setSingleStep(0.010);
    m_saltConcSpin->setDecimals(3);
    m_saltConcSpin->setValue(0.050);
    m_saltConcSpin->setSuffix("  M");
    saltForm->addRow("Concentration:", m_saltConcSpin);
    m_ionStack->addWidget(saltPage);

    // Page 1
    auto *manualPage = new QWidget;
    auto *manualForm = new QFormLayout(manualPage);
    manualForm->setContentsMargins(0, 4, 0, 0);
    m_posIonsSpin = new QSpinBox;
    m_posIonsSpin->setRange(0, 500);
    m_posIonsSpin->setValue(0);
    m_negIonsSpin = new QSpinBox;
    m_negIonsSpin->setRange(0, 500);
    m_negIonsSpin->setValue(0);
    manualForm->addRow("Positive ions:", m_posIonsSpin);
    manualForm->addRow("Negative ions:", m_negIonsSpin);
    m_ionStack->addWidget(manualPage);

    ionVbox->addWidget(m_ionStack);

    // Mg2+ toggle applies to both modes
    m_magnesiumCheck = new QCheckBox("Use Mg²⁺ as positive ion  (instead of Na⁺)", this);
    m_magnesiumCheck->setChecked(false);
    ionVbox->addWidget(m_magnesiumCheck);

    root->addWidget(ionGroup); // add to root layout, not the form

    // Switch pages when radio changes
    connect(m_ionModeGroup, &QButtonGroup::idClicked,
            m_ionStack, &QStackedWidget::setCurrentIndex);

    // Simulation time
    m_simTimeSpin = new QDoubleSpinBox(this);
    m_simTimeSpin->setRange(0.1, 1000.0);
    m_simTimeSpin->setSingleStep(10.0);
    m_simTimeSpin->setDecimals(1);
    m_simTimeSpin->setValue(100.0);
    m_simTimeSpin->setSuffix("  ns");
    form->addRow("Simulation time:", m_simTimeSpin);

    root->addWidget(simGroup);

    // ── Run button ───────────────────────────────────────────────────────
    auto *runRow = new QHBoxLayout;
    runRow->addStretch(1);

    m_stopButton = new QPushButton("■  Stop", this);
    m_stopButton->setObjectName("stopButton");
    m_stopButton->setEnabled(false);
    connect(m_stopButton, &QPushButton::clicked, this, &MainWindow::onStopClicked);
    runRow->addWidget(m_stopButton);

    m_runButton = new QPushButton("▶  Run simulation", this);
    m_runButton->setObjectName("runButton");
    connect(m_runButton, &QPushButton::clicked, this, &MainWindow::onRunClicked);
    runRow->addWidget(m_runButton);

    root->addLayout(runRow);

    // ── Progress bar ─────────────────────────────────────────────────────
    m_progressBar = new QProgressBar(this);
    m_progressBar->setRange(0, 100);
    m_progressBar->setValue(0);
    m_progressBar->setTextVisible(true);
    m_progressBar->setFormat("Production MD — %p%");
    m_progressBar->setVisible(false); // hidden until simulation starts
    root->addWidget(m_progressBar);

    // ── Log output ───────────────────────────────────────────────────────
    auto *logGroup = new QGroupBox("Output log", this);
    auto *logLayout = new QVBoxLayout(logGroup);
    logLayout->setContentsMargins(8, 12, 8, 8);

    m_logOutput = new QTextEdit(this);
    m_logOutput->setReadOnly(true);
    m_logOutput->setMinimumHeight(160);
    logLayout->addWidget(m_logOutput);
    root->addWidget(logGroup, 1);
}

// ── slots ──────────────────────────────────────────────────────────────────
void MainWindow::onBrowseScript()
{
    QSettings settings;
    const QString lastDir = settings.value("lastScriptDir", QDir::homePath()).toString();

    const QString path = QFileDialog::getOpenFileName(
        this, "Select chikaterasu.sh", lastDir,
        "Shell scripts (*.sh);;All files (*)");

    if (!path.isEmpty())
    {
        m_scriptPathEdit->setText(path);
        settings.setValue("lastScriptDir", QFileInfo(path).absolutePath());
    }
}

void MainWindow::onRunClicked()
{
    if (scriptPath().isEmpty())
    {
        QMessageBox::warning(this, "No script", "Please select the chikaterasu.sh script first.");
        return;
    }

    m_progressBar->setValue(0);
    m_progressBar->setVisible(true);

    m_totalSteps = static_cast<int>(m_simTimeSpin->value() * 500000);
    buildAndWriteConfig();

    m_logOutput->clear();
    m_logOutput->append("<span style='color:#4299e1'>Starting simulation…</span>");

    m_runButton->setEnabled(false);
    m_runButton->setText("Running…");
    m_stopButton->setEnabled(true);

    // Patch env variables and call the script
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("CHIKA_GUI", "1");
    env.insert("CHIKA_PROTEIN", m_proteinNameEdit->text().trimmed());
    env.insert("CHIKA_TEMP", QString::number(m_temperatureSpin->value(), 'f', 1));
    env.insert("CHIKA_SALTCONC", QString::number(m_saltConcSpin->value(), 'f', 3));
    env.insert("CHIKA_SIMTIME", QString::number(m_simTimeSpin->value(), 'f', 1));
    m_process->setProcessEnvironment(env);
    m_process->setWorkingDirectory(QFileInfo(scriptPath()).absolutePath());

    m_process->setProgram("setsid");
    m_process->setArguments({"bash", "-c",
                             "source /usr/local/gromacs/bin/GMXRC && bash \"" + scriptPath() + "\""});
    m_process->start();
}

void MainWindow::onProcessOutput()
{
    const QString out = m_process->readAllStandardOutput();
    const QString err = m_process->readAllStandardError();

    if (!out.isEmpty())
        m_logOutput->append(out.trimmed());

    if (!err.isEmpty())
    {
        static QRegularExpression reStep(R"(^step\s+(\d+),)");
        const QStringList lines = err.split(QRegularExpression("[\r\n]"), Qt::SkipEmptyParts);

        for (const QString &line : lines)
        {
            QRegularExpressionMatch match = reStep.match(line.trimmed());
            if (match.hasMatch() && m_totalSteps > 0)
            {
                const int step = match.captured(1).toInt();
                if (step > 0)
                    m_progressBar->setValue(step * 100 / m_totalSteps);
            }
            else
            {
                const QString trimmed = line.trimmed();
                if (!trimmed.isEmpty())
                    m_logOutput->append(
                        "<span style='color:#fc8181'>" + trimmed.toHtmlEscaped() + "</span>");
            }
        }
    }

    m_logOutput->verticalScrollBar()->setValue(
        m_logOutput->verticalScrollBar()->maximum());
}

void MainWindow::onProcessFinished(int exitCode)
{
    m_progressBar->setValue(100);
    m_runButton->setEnabled(true);
    m_runButton->setText("▶  Run simulation");
    m_stopButton->setEnabled(false);
    if (exitCode == 0)
        m_logOutput->append("<span style='color:#68d391'> Finished successfully.</span>");
    else
        m_logOutput->append(
            QString("<span style='color:#fc8181'> Exited with code %1.</span>").arg(exitCode));
}

void MainWindow::onStopClicked()
{
    if (m_process->state() == QProcess::NotRunning)
        return;

    const qint64 pid = m_process->processId();
    if (pid > 0)
    {
        // Kill everything in this process's session
        QProcess::execute("pkill", {"-TERM", "-s", QString::number(pid)});
        // Fallback: also kill the process group directly
        QProcess::execute("kill", {"-TERM", "-" + QString::number(pid)});
    }

    m_logOutput->append("<span style='color:#fc8181'>⏹ Stopped by user.</span>");
    m_process->kill(); // ensure QProcess itself is also reaped
}

// ── helpers ────────────────────────────────────────────────────────────────
QString MainWindow::scriptPath() const
{
    return m_scriptPathEdit->text().trimmed();
}

/**
 * Write a small config file next to the script so the bash script can
 * source it instead of having hard-coded values at the top.
 * The bash script just needs:
 *   source "$(dirname "$0")/chika_gui.conf"
 */
void MainWindow::buildAndWriteConfig()
{
    const QString dir = QFileInfo(scriptPath()).absolutePath();
    const QString conf = dir + "/chika_gui.conf";

    QFile f(conf);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream ts(&f);
    ts << "# Auto-generated by ChikaterasuGUI — do not edit manually\n";
    ts << "protein_name=\"" << m_proteinNameEdit->text().trimmed() << "\"\n";
    ts << "ref_t=" << QString::number(m_temperatureSpin->value(), 'f', 1) << "\n";
    ts << "sim_time_ns=" << QString::number(m_simTimeSpin->value(), 'f', 1) << "\n";
    ts << "debug_level=" << m_debugLevelCombo->currentData().toInt() << "\n";
    // ts << "his_manual=" << (m_hisManualCheck->isChecked() ? "true" : "false") << "\n";
    ts << "water=" << m_waterModelCombo->currentData().toString() << "\n";
    ts << "nruns=" << m_nrunsSpin->value() << "\n";
    ts << "disulfide=" << (m_disulfideCheck->isChecked() ? "true" : "false") << "\n";
    ts << "forcefield=" << m_forceFieldCombo->currentData().toString() << "\n";

    const bool useSalt = m_ionModeSalt->isChecked();
    ts << "specify_salt_concentration=" << (useSalt ? "true" : "false") << "\n";
    ts << "salt_concentration=" << QString::number(m_saltConcSpin->value(), 'f', 3) << "\n";
    ts << "pos_ions=" << m_posIonsSpin->value() << "\n";
    ts << "neg_ions=" << m_negIonsSpin->value() << "\n";
    ts << "magnesium=" << (m_magnesiumCheck->isChecked() ? "true" : "false") << "\n";

    const bool boxManual = m_boxManualCheck->isChecked();
    ts << "box_manual=" << (boxManual ? "true" : "false") << "\n";
    ts << "box_dim=\"    " << m_boxDimEdit->text().trimmed() << " \"\n";
    ts << "cell_shape=\"" << m_cellShapeCombo->currentData().toString() << "\"\n";

    ts << "shear_enabled=" << (m_shearCheck->isChecked() ? "true" : "false") << "\n";
    ts << "shear_rate=" << QString::number(m_shearRateSpin->value(), 'f', 4) << "\n";
}
